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

# These are simply wrappers around strings, that remove whitespace
for T in (:Sample, :Clade)
    @eval begin
        struct $T
            name::String

            function $T(s::Union{String, SubString{String}})
                str = strip(s)
                isempty(str) && error("String must be nonempty!")
                if str == s
                    new(s isa String ? s : String(s))
                else
                    new(String(str))
                end
            end
        end

        function Base.tryparse(::Type{$T}, s::AbstractString)
            str = strip(s)
            isempty(str) ? nothing : $T(str)
        end

        Base.parse(::Type{$T}, s::AbstractString) = $T(s)

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
split_clade(s) = parseout_suffix(Clade, s, '_')

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
include("blast.jl")

# This is type piracy, but hopefully it's OK!
StructTypes.StructType(::Type{<:BioSequence}) = StructTypes.StringType()
StructTypes.StructType(::Type{<:AbstractUnitRange}) = StructTypes.Struct()

"""
    store_references(dst::Union{IO, AbstractString}, ref_itr)

Serialize to file at path or IO `dst` an iterable of `Reference}.
"""
function store_references end

function store_references(io::IO, reference_itr)
    JSON3.write(io, vec(collect(reference_itr))::Vector{Reference})
end

function store_references(path::AbstractString, reference_itr)
    open(path, "w") do io
        store_references(io, reference_itr)
    end
end

"""
    load_references(file::AbstractString)

Loads a JSON file created by `store_references`. If the version number of the
Influenza package used to serialize the JSON is not compatible with the version
number of `Influenza` used to load, throw an error.
"""
function load_references(file::AbstractString)
    open(file) do io
        JSON3.read(io, Vector{Reference})
    end
end

export TERMINAL_INFLUENZA_5,
    TERMINAL_INFLUENZA_3,

    Sample,
    Clade,
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
    parseout_suffix,
    split_segment,
    split_protein,
    split_clade,

    # Exports from InfluenzaCore
    Segment, Segments, SeroType, Proteins, Protein, source

end # module
