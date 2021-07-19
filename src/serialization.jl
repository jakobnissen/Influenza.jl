# This file contains to code write References to .jls files and .fna files.

"""
    store_references(dir::AbstractString, ref_itr, fasta::Bool=false)

For each segment in the iterator of `Segment`, `ref_itr`, serialize first the
current Influenza version as `VersionNumber`then a `Vector{Reference}`.
If `fasta`, also write a FASTA file to `\$segment.fna`
"""
function store_references(
    dir::AbstractString,
    reference_itr,
    fasta::Bool=false,
)
    isdir(dir) || error("No such directory: \"$dir\"")

    # Collect by segment
    bysegment = Dict{Segment, Vector{Reference}}()
    for ref::Reference in reference_itr
        push!(get!(Vector{Reference}, bysegment, ref.segment), ref)
    end

    # Serialize each segment
    for (segment, refs) in bysegment
        open(joinpath(dir, "$segment.jls")) do io
            serialize(io, VERSION)
            serialize(io, refs)
        end

        if fasta
            open(FASTA.Writer, joinpath(dir, "$segment.fna")) do writer
                for ref in refs
                    write(writer, ref.seq)
                end
            end
        end
    end
end

function is_breaking(a::VersionNumber, b::VersionNumber)
    return a.major != b.major || (iszero(a.major) && (a.minor != b.minor))
end

"""
    load_references(file::AbstractString)

Loads a `.jls` file created by `store_references`. If the version number of the
`jls` file is not compatible with the version number of `Influenza` used to load,
throw an error.
"""
function load_references(file::AbstractString)
    open(file) do io
        v = deserialize(io)
        if is_breaking(v, INFLUENZA_VERSION)
            error(
                "File \"$file\" was serialized with Influenza version $v, ",
                "cannot load with current Influenza version $INFLUENZA_VERSION."
            )
        end
        deserialize(io)
    end
end


