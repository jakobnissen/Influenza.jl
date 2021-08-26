# This file contains to code write References to .jls files and .fna files.
"""
    store_references(basename::AbstractString, ref_itr, fasta::Bool=false)

Serialize to file `basename * ".jls"` first the current Influenza version as
`VersionNumber`, then a `Vector{Reference}`. If `fasta`, also write a FASTA file
to `basename * ".fna"`.
"""
function store_references(
    basename::AbstractString,
    reference_itr,
    fasta::Bool=false,
)
    # Check existence of parent directory
    parent_dir = dirname(basename)
    if !isempty(parent_dir) && !isdir(parent_dir)
        error("Parent directory not found: \"$parent_dir\"")
    end

    # Serialize references
    ref_vec::Vector{Influenza.Reference} = collect(reference_itr)
    open(basename * ".jls", "w") do io
        serialize(io, INFLUENZA_VERSION)
        serialize(io, ref_vec)
    end

    # Serialize as FASTA format
    if fasta
        open(FASTA.Writer, basename * ".fna") do writer
            for ref in ref_vec
                header = ref.name * '_' * string(ref.segment)
                write(writer, FASTA.Record(header, ref.seq))
            end
        end
    end
end

function is_breaking(a::VersionNumber, b::VersionNumber)
    return a.major != b.major || (iszero(a.major) && (a.minor != b.minor))
end

"Strip <SEP><SEGMENT> off the end of a string. Errors if not present."
function strip_trailing_segment(s::Union{String, SubString{String}}, sep::UInt8=UInt8('_'))
    p = let
        _p = findlast(isequal(sep), codeunits(s))
        _p === nothing ? error("Sep '", Char(sep), "' not found in header \"", s, '"') : _p
    end
    lastpart = view(s, p+1:lastindex(s))
    if tryparse(Segment, lastpart) === nothing
        error("Cannot parse as segment: \"", lastpart, '"')
    end
    return view(s, 1:prevind(s, p))
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
        deserialize(io)::Vector{Reference}
    end
end


