# Code to analyze new sequences using BLAST and serialized References.
tovec(x) = x isa Vector ? x : vec(collect(x))

"""
    annotate(itr, jls_path, fna_path) -> Vector{Option{AlignedAssembly}}

Annotate `itr`, a vector of `FASTA.Record` or iterator of `Assembly`, yielding
a vector of `Option{AlignedAssembly}`, with error values where no reference could
be matched. `jls_path` and `fna_path` must be paths to serialization and FASTA
file produced by `Influenza.store_references.`

# Extended help
This function works by using running `blastn` in a temporary directory and
parsing the output. Hence, `blastn` must be available on Julia's `PATH`. It is
fairly slow for hundreds of sequences, using >100 ms per seq.
"""
function annotate end

function annotate(
    itr,
    jls_path::AbstractString,
    fna_path::AbstractString,
)
    # Check presence of files
    for path in (jls_path, fna_path)
        isfile(path) || error("File not found: \"$path\"")
    end

    # Collect to vector of assembly
    assembly_vec::Vector{Assembly} = tovec(itr)

    # Write FASTA to temp file. We name the records by index, e.g. 1, 2, 3...
    tmp_dir = tempdir()
    (tmp_path, tmp_file) = mktemp(tmp_dir)
    writer = FASTA.Writer(tmp_file)
    for (i, assembly::Assembly) in pairs(assembly_vec)
        write(writer, FASTA.Record(string(i), assembly.seq))
    end
    close(writer)

    # Run blastn
    blast_path = tempname(tmp_dir)
    command = `blastn -query $tmp_path -subject $fna_path -outfmt "6 qacc sacc bitscore qlen length pident"`
    run(pipeline(command, stdout=blast_path))

    # Parse blastn output to Dict{String, Vector{Int}} of subject => [n_assembly ... ]
    best_hits::Dict{String, Vector{Int}} = let
        maybe_best_hits::Dict{Int, Option{String}} = open(blast_path) do io
            parse_blastout(io, 0.9)
        end
        res = Dict{String, Vector{Int}}()
        for i in eachindex(assembly_vec)
            subject = @unwrap_or get(maybe_best_hits, i, none(String)) continue
            push!(get!(Vector{Int}, res, subject), i)
        end
        res
    end

    # Load references
    loaded_refs = load_references(jls_path)
    references = fill(none(Reference), length(assembly_vec))
    for ref in loaded_refs
        vecnums = get(best_hits, ref.name, nothing)
        if vecnums === nothing
            continue
        else
            for i in vecnums
                references[i] = some(ref)
            end
        end
    end

    map(zip(assembly_vec, references)) do (assembly, maybe_ref)
        and_then(ref -> AlignedAssembly(assembly, ref), AlignedAssembly, maybe_ref)
    end
end

function annotate(
    recs::Vector{FASTA.Record},
    jls_path::AbstractString,
    fna_path::AbstractString,
)
    assemblies = map(recs) do rec
        Assembly(rec, nothing)
    end
    validate(assemblies, serialize_dir)
end

function parse_blastout(io::IO, lenratio::Real)::Dict{Int, Option{String}}
    # Parse to vector of hits
    hits = eachline(io) |> imap(strip) |> ifilter(!isempty) |> imap() do line
        fields = split(line, '\t')
        @assert length(fields) == 6
        query = parse(Int, fields[1])
        subject = String(fields[2])
        bitscore = parse(Float64, fields[3])
        qlen = parse(UInt, fields[4])
        len = parse(UInt, fields[5])
        ident = parse(Float64, fields[6]) / 100
        (; query, subject, bitscore, qlen, len, ident)
    end |> collect

    # Group by query sequence, sort by bitscore
    byquery = Dict{Int, Vector{eltype(hits)}}()
    for hit in hits
        push!(get!(Vector{eltype(hits)}, byquery, hit.query), hit)
    end

    foreach(values(byquery)) do hits
        sort!(hits, by=x -> x.bitscore, rev=true)
    end

    # Get best hits
    result = Dict{Int, Option{String}}()
    for (query, hits) in byquery
        for hit in hits
            if hit.len / hit.qlen ≥ lenratio && hit.ident ≥ 0.9
                result[hit.query] = some(hit.subject)
                break
            end
        haskey(result, query) || (result[query] = none(String))
        end
    end
    return result
end
