"""
    InfluenzaCore

This module defines only *core types* for the other Influenza packages.
It is *not* meant to be external, application-faced, nor is it meant to include
functions and algorithms for interactive use.

Its only function is to be imported by scripts and such in order to get access
to the most basic types.
"""
module InfluenzaCore

baremodule Segments
import Base: @enum, @doc
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
Base.tryparse(::Type{Segment}, s::AbstractString) = get(STRING_SEGMENT_DICT, strip(s), nothing)

baremodule Proteins
import Base: @enum, @doc, Dict, =>
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

Base.tryparse(::Type{Protein}, s::AbstractString) = get(Proteins._STR_PROTEINVARIANT, strip(s), nothing)

const PROTEIN_TO_SEGMENT = Tuple(UInt8[0,1,1,1,2,2,3,4,5,5,6,6,6,6,7,7,7])
@assert length(PROTEIN_TO_SEGMENT) == length(instances(Protein))
@assert maximum(PROTEIN_TO_SEGMENT) == length(instances(Segment)) - 1

function source(x::Protein)
    index = reinterpret(UInt8, x) + 0x01
    integer = @inbounds PROTEIN_TO_SEGMENT[index]
    return reinterpret(Segment, integer)
end

"""
SeroType

Influenza A are separated into serotypes, based on their HA and NA proteins.
The H and Ns are indicated by a natural number (1, 2, 3...), or may optionally
be `nothing` (missing types conventionally written as 0).
"""
struct SeroType
    h::Union{Nothing, UInt8}
    n::Union{Nothing, UInt8}

    function SeroType(h_, n_)
        h = convert(Union{Nothing, UInt8}, h_)
        n = convert(Union{Nothing, UInt8}, n_)
        if (h !== nothing && iszero(h)) || (n !== nothing && iszero(n))
            throw(DomainError(0x00, "Serotypes must be nonzero"))
        end
        new(h, n)
    end
end

function Base.show(io::IO, x::SeroType)
    dummy(x) = isnothing(x) ? '0' : x
    print(io, 'H', dummy(x.h), 'N', dummy(x.n))
    return nothing
end

function Base.tryparse(::Type{SeroType}, x::AbstractString)
    m = match(r"^H(\d+)N(\d+)$", x)
    m === nothing && return nothing
    H = parse(UInt8, m[1]::SubString{String})
    H = iszero(H) ? nothing : H
    N = parse(UInt8, m[2]::SubString{String})
    N = iszero(N) ? nothing : N
    return SeroType(H, N)
end

export Segment, Segments, SeroType, Proteins, Protein, source

end # module
