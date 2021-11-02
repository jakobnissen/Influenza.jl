"""
This module contains generally useful functionality related to influenza processing. It is imported by other packages.

It is distinct from `InfluenzaCore` in that it is a more heavyweight dependency, and distinct from other packages in that this package does not serve any particular purpose, but as a swiss knife tool for Influenza analysis.
"""
module Influenza

using ErrorTypes: Option, some, none, is_error, unwrap, unwrap_or, @unwrap_or, and_then
using FASTX: FASTA
using InfluenzaCore: Segment, Protein, Segments, SeroType, Proteins, source
using StructTypes: StructTypes # for JSON (de)serialization
using JSON3: JSON3
using BioSequences: BioSequences, LongDNASeq, LongAminoAcidSeq, NucleotideSeq,
    BioSequence, DNACodon, @mer_str, isgap, @aa_str, @biore_str, NucleicAcidAlphabet,
    DNA, DNA_Gap
using BioAlignments: BioAlignments
using Printf: Printf, @sprintf
using BlastParse: BlastParse

const BA = BioAlignments

imap(f) = x -> Iterators.map(f, x)
ifilter(f) = x -> Iterators.filter(f, x)

"Conserved sequence at influenza segment 5' end"
const TERMINAL_INFLUENZA_5 = mer"AGCAAAAGCAGG"dna

"Conserved sequence at influenza segment 3' end"
const TERMINAL_INFLUENZA_3 = mer"CTTGTTTCTCCT"dna

const INFLUENZA_VERSION = let
    p = pathof(Influenza)::String
    project = joinpath(dirname(dirname(p)), "Project.toml")
    toml = read(project, String)
    m = match(r"(*ANYCRLF)^version\s*=\s\"(.*)\"$"m, toml)::RegexMatch
    VersionNumber(m[1]::SubString{String})
end

# These are simply wrappers around strings, that remove whitespace
for T in (:Sample,)
    @eval begin
        struct $T
            name::String

            function $T(s::Union{String, SubString{String}})
                str = strip(s)
                if str == s
                    new(s isa String ? s : String(s))
                else
                    new(String(str))
                end
            end
        end

        function $T(s::AbstractString)
            $T(convert(String, s))
        end

        Base.nameof(x::$T) = x.name
        Base.hash(x::$T, h::UInt) = hash(nameof(x), hash($T, h))
        Base.:(==)(x::$T, y::$T) = nameof(x) == nameof(y)
        Base.print(io::IO, x::$T) = print(io, nameof(x))
        Base.isless(x::$T, y::$T) = isless(nameof(x), nameof(y))
    end
end

"""
    try_parseout_suffix(::Type{T}, s::AbstractString, sep::Char)

Given a string `s` of format `X * sep * Y`, return (`X`, parse(T, Y))`, or `nothing`
if the input string cannot be parsed as that format.

# Example
```julia
julia> try_parseout_suffix(Segment, "Denmark_PB1", '_')
("Denmark", InfluenzaCore.Segments.PB1)

julia> try_parseout_suffix(Segment, "Denmark_WWx", '_') === nothing
true

julia> try_parseout_suffix(Segment, "xxx", '_') === nothing
true
```
"""
function try_parseout_suffix(::Type{T}, s::AbstractString, sep::Char) where T
    p = findlast(==(sep), s)
    p === nothing && return nothing
    val = tryparse(T, SubString(s, nextind(s, p):lastindex(s)))
    val === nothing && return nothing
    return (SubString(s, 1, prevind(s, p)), val)
end

"Like `try_parseout_suffix`, but errors instead of returning `nothing`"
function parseout_suffix(::Type{T}, s::AbstractString, sep::Char) where T
    pair = try_parseout_suffix(T, s, sep)
    if pair === nothing
        error("Expected String * $sep * $T, got \"", s, '\"')
    else
        pair
    end
end

split_segment(s) = parseout_suffix(Segment, s, '_')
split_protein(s) = parseout_suffix(Protein, s, '_')

"""
    groupby(f, itr)

Return a `Dict` where the keys are the unique values of `f` applied to each element
of `itr`
"""
function groupby(f, itr)
    tupelems(::Vector{Tuple{A, B}}) where {A, B} = (A, B)
    v = [(f(i), i) for i in itr]
    K, V = tupelems(v)
    d = Dict{K, Vector{V}}()
    for (k, val) in v
        push!(get!(valtype(d), d, k), val)
    end
    d
end

"Singleton struct used for three-valued logic, in e.g. `Union{Bool, Maybe}`"
struct Maybe end

include("errors.jl")
include("alignment.jl")
include("assembly.jl")
include("serialization.jl")
include("blast.jl")

# This is type piracy, but hopefully it's OK!
StructTypes.StructType(::Type{<:BioSequence}) = StructTypes.StringType()
StructTypes.StructType(::Type{<:AbstractUnitRange}) = StructTypes.Struct()

export TERMINAL_INFLUENZA_5,
    TERMINAL_INFLUENZA_3,

    Sample,
    Assembly,
    AssemblyProtein,
    Reference,
    AlignedAssembly,

    groupby,

    translate_reference,
    translate_proteins,
    alignment_identity,
    ha0_cleavage,
    store_references,
    load_references,
    annotate,

    try_parseout_suffix,
    parseout_suffix

    # Exports from InfluenzaCore
    Segment, Segments, SeroType, Proteins, Protein, source

end # module
