# Code to analyze new sequences using BLAST and serialized References.
tovec(x) = x isa Vector ? x : vec(collect(x))

"""
    annotate(itr, json_path, subject, mmpath) -> Vector{Option{AlignedAssembly}}

Annotate `itr`, a vector of `FASTA.Record` or iterator of `Assembly`, yielding
a vector of `Option{AlignedAssembly}`, with error values where no reference could
be matched. `json_path` must be path to serialization and FASTA
file produced by `Influenza.store_references.` The `subject` can be a FASTA file
of the same sequences as the JSON, or an mmseqs database. `mmseqs` can be `nothing`,
or the path of the mmseqs `tmp` directory, for speed.
It is recommended that then references broadly represent various clades,
so a relatively close hit can be found.

# Extended help
This function works by using running `blastn` in a temporary directory and
parsing the output. Hence, `blastn` must be available on Julia's `PATH`. It is
fairly slow for hundreds of sequences, using >100 ms per seq.
"""
function annotate end

# Convert records to vector of assembly, then dispatch
function annotate(recs::Vector{FASTA.Record}, args...)
    assemblies = map(recs) do rec
        Assembly(rec, nothing)
    end
    annotate(assemblies, args...)
end

# Convert iterable to vector of assembly, then dispatch
function annotate(itr, args...)
    assemblies = tovec(itr)::Vector{Assembly}
    annotate(assemblies, args...)
end

# Check paths
function annotate(
    itr::Vector{Assembly},
    json_path::AbstractString,
    fna_path::AbstractString,
    mm_tmpdir::Union{Nothing, AbstractString}=nothing,
)
    # Check presence of files
    for path in (json_path, fna_path, mm_tmpdir)
        if !isnothing(path) && !ispath(path)
            error("File not found: \"$path\"")
        end
    end

    tmp_dir = mktempdir()
    mm_tmpdir = if mm_tmpdir === nothing
        mktempdir(tmp_dir)
    else
        mm_tmpdir
    end

    _annotate(itr, json_path, fna_path, tmp_dir, mm_tmpdir)
end

function _annotate(
    assembly_vec::Vector{Assembly},
    json_path::AbstractString,
    fna_path::AbstractString,
    tmp_dir::AbstractString,
    mm_tmpdir::AbstractString,
)
    # Check presence of mmseqs
    run(pipeline(`which mmseqs`, stdout=devnull))

    # Write FASTA to temp file. We name the records by index, e.g. 1, 2, 3...
    (tmp_path, tmp_file) = mktemp(tmp_dir)
    writer = FASTA.Writer(tmp_file)
    for (i, assembly::Assembly) in pairs(assembly_vec)
        write(writer, FASTA.Record(string(i), assembly.seq))
    end
    close(writer)

    # Run blastn
    blast_path = tempname(tmp_dir)
    mmseqs_search(tmp_path, fna_path, blast_path, mm_tmpdir, tempname(tmp_dir))

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
    loaded_refs = load_references(json_path)
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

    result = fill(none(AlignedAssembly), length(references))
    Threads.@threads for i in eachindex(result)
        asm = assembly_vec[i]
        mref = references[i]
        v = and_then(ref -> AlignedAssembly(asm, ref), AlignedAssembly, mref)
        result[i] = v
    end
    return result
end

function mmseqs_search(
    query::AbstractString,
    target::AbstractString,
    result::AbstractString,
    tmpdir::AbstractString,
    stdout::AbstractString
)
    # Search type 3 is nucleotide, the rest means "fast but insensitive".
    cmd = `mmseqs easy-search $query $target $result $tmpdir
    --search-type 3 --exact-kmer-matching 1 -s 1
    --format-output "query,target,bits,qlen,alnlen,pident"`
    run(pipeline(cmd, stdout=stdout))
end

@eval $(BlastParse.gen_blastparse_code(
    (:qacc, :sacc, :bitscore, :qlen, :length, :pident),
    :parse_blast_io
))

function parse_blastout(io::IO, lenratio::Real)::Dict{Int, Option{String}}
    # Parse to vector of hits
    hits = parse_blast_io(io)

    # Group by query sequence, sort by bitscore
    byquery = Dict{Int, Vector{eltype(hits)}}()
    for hit in hits
        push!(get!(valtype(byquery), byquery, parse(Int, hit.qacc)), hit)
    end

    foreach(values(byquery)) do hits
        sort!(hits, by=x -> x.bitscore, rev=true)
    end

    # Get best hits
    result = Dict{Int, Option{String}}()
    for (query, hits) in byquery
        for hit in hits
            if hit.length / hit.qlen ≥ lenratio && hit.pident ≥ 0.8
                # The result has format NAME_SEGMENT, but the refs in the json file
                # does not have the trailing segment, so we strip it off here.
                result[parse(Int, hit.qacc)] = some(String(hit.sacc))
                break
            end
        haskey(result, query) || (result[query] = none(String))
        end
    end
    return result
end
