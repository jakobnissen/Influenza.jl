"""
    InfluenzaCore

This module defines only *core types* for the other Influenza packages.
It is *not* meant to be external, application-faced, nor is it meant to include
functions and algorithms for interactive use.

Its only function is to be imported by scripts and such in order to get access
to the most basic types.
"""
module InfluenzaCore

using ErrorTypes
using SumTypes

module Segments
"""
    Segment

An `Enum` representing an Influenza genomic segment. It is ordered from longest (PB2)
to shortest (NS), as per canonical ordering.
"""
@enum Segment::UInt8 PB2 PB1 PA HA NP NA MP NS
export Segment
end
using .Segments

const STRING_SEGMENT_DICT = Dict(string(s)=>s for s in instances(Segment))
function Base.parse(::Type{Segment}, s::AbstractString)::Option{Segment}
    y = get(STRING_SEGMENT_DICT, s, nothing)
    y === nothing && return none
    some(y)
end

module Proteins
"""
    Protein

An `Enum` representing an influenza protein. Not all virions have all proteins, some are
auxiliary or specific to some species of influenza.
"""
@enum Protein::UInt8 begin
    PB2
    PB1
    N40 # very rarely seen
    PB1F2 # inf A only
    PA
    PAX # inf A only
    HA
    NP
    NA
    NB # inf B only
    M1
    M2 # inf A only
    BM2 # inf B only
    M42 # uncommon
    NS1
    NEP
    NS3 # uncommon
end

const _STR_PROTEINVARIANT = Dict(
    "PB2" => PB2,
    "PB1" => PB1,
    "N40" => N40,
    "PB1-F2" => PB1F2,
    "PB1F2" => PB1F2,
    "PB1-FA" => PB1F2,
    "PA" => PA,
    "PA-X" => PAX,
    "PAX" => PAX,
    "HA" => HA,
    "NP" => NP,
    "NA" => NA,
    "NB" => NB,
    "M1" => M1,
    "M2" => M2,
    "BM2" => BM2,
    "M42" => M42,
    "NS1" => NS1,
    "NS2" => NEP,
    "NEP" => NEP, # yeah, it has two names
    "NS3" => NS3,
)

export Protein
end
using .Proteins

function Base.parse(::Type{Protein}, s::AbstractString)::Option{Protein}
    y = get(Proteins._STR_PROTEINVARIANT, s, nothing)
    y === nothing ? none : some(y)
end

const PROTEIN_TO_SEGMENT = Tuple(UInt8[0,1,1,1,2,2,3,4,5,5,6,6,6,6,7,7,7])
@assert length(PROTEIN_TO_SEGMENT) == length(instances(Protein))
@assert maximum(PROTEIN_TO_SEGMENT) == length(instances(Segment)) - 1

function source(x::Protein)
    index = reinterpret(UInt8, x) + 0x01
    integer = @inbounds PROTEIN_TO_SEGMENT[index]
    return reinterpret(Segment, integer)
end

module SubTypes
using SumTypes

@sum_type SubType begin
    Yamagata()
    Victoria()
    InfluenzaA(::UInt8, ::UInt8)
end

InfluenzaA(H::Integer, N::Integer) = InfluenzaA(UInt8(H), UInt8(N))
export SubType
end
using .SubTypes

Base.show(io::IO, x::SubTypes.InfluenzaA) = print(io, "H$(x._1)N$(x._2)")

function Base.parse(::Type{SubType}, s::AbstractString)::Option{SubType}
    s == "Victoria" && return some(SubTypes.Victoria())
    s == "Yamagata" && return some(SubTypes.Yamagata())
    m = match(r"^H(\d+)N(\d+)$", s)
    m === nothing && return none
    H = parse(UInt8, m[1]::SubString{String})
    N = parse(UInt8, m[2]::SubString{String})
    some(SubTypes.InfluenzaA(H, N))
end

export Segment, Segments, SubType, SubTypes, Proteins, Protein, source

end # module
