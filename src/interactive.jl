BLASTN = "/Users/jakobnissen/miniconda3/bin/blastn"
REFDIR = "/Users/jakobnissen/Documents/ssi/projects/flupipe/ref/seqs"

"""
    manual_check(::Segment, ::LongDNASeq; refdir=[preset])

Checks a segment sequence by BLASTing it to a collection of references.
Returns a (ReferenceAssembly, Vector{String}) result, with the strings being
the lines would be in a flupipe report.

Note that this approach is fairly inefficient, and should probably not be used in
a loop unless you want your CPU to keep you warm in the winter.
"""
function manual_check(segment::Segment, seq::LongDNASeq, refdir::String=REFDIR)
    # Write to file
    filename = tempname()
    open(filename, "w") do file
        println(file, ">seq\n", seq)
    end

    # BLAST to ref
    outfile = tempname()
    subject = joinpath(refdir, "$segment.fna")
    command = `$BLASTN -query $filename -subject $subject -outfmt "6 sacc bitscore length pident"`
    run(pipeline(command, stdout=outfile))


    best_hit = open(outfile) do file
        x = parse_blastout(file, trunc(Int, 0.9 * length(seq)))
        expect(x, "No good hits found")
    end

    # Extract best hit, get that hit from jls file
    reference = load_references(segment, refdir, Set([best_hit]))[best_hit]
    
    # Create ReferenceAssembly
    assembly = Assembly(segment, falses(length(seq)), seq, best_hit)
    refasm = ReferenceAssembly(reference, assembly)
    return (refasm, segment_report_lines(segment, some(refasm), Dict{String, Any}()))
end

function parse_blastout(io::IO, minlength::Integer)::Option{String}
    fields = Vector{SubString{String}}(undef, 4)
    hits = eachline(io) |> Map(strip) |> Filter(!isempty) |> Map() do line
        split!(fields, line, UInt8('\t'))
        bitscore = parse(UInt, fields[2])
        len = parse(UInt, fields[3])
        ident = parse(Float64, fields[4]) / 100
        (String(first(fields)), bitscore, len, ident)
    end |> collect
    sort!(hits, by=x -> x[2], rev=true)

    for (name, bitscore, len, ident) in hits
        if len ≥ minlength && ident ≥ 0.9
            return some(name)
        end
    end
    return none
end