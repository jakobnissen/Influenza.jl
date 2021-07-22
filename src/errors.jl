"Struct representing an indel compared to the reference"
struct Indel
    # range: Position in ref (del) or asm (ins) of bases affected
    range::UnitRange{UInt32}
    # position: seq aligns between pos and pos+1 in the other seq
    position::UInt32
    is_deletion::Bool

    function Indel(range, pos, isdel)
        rng = convert(UnitRange{UInt32}, range)
        isempty(rng) && throw(ArgumentError("Cannot have zero-length indel"))
        new(rng, convert(UInt32, pos), convert(Bool, isdel))
    end
end

Base.length(x::Indel) = length(x.range)

function indel_message(x::Indel)
    rangestring = string(first(x.range)) * '-' * string(last(x.range))
    posstring = string(x.position) * '/' * string(x.position + 1)
    if x.is_deletion
        "Deletion of ref pos " * rangestring * " b/w pos " * posstring
    else
        "Insertion of bases " * rangestring * " b/w ref pos " * posstring
    end
end

abstract type InfluenzaError end
abstract type SegmentError <: InfluenzaError end
abstract type ProteinError <: InfluenzaError end

"Protein or DNA sequence has too low identity compared to the reference"
struct ErrorLowIdentity <: InfluenzaError
    identity::Float32
end

function Base.print(io::IO, x::ErrorLowIdentity)
    percent = round(x.identity * 100, digits=1)
    print(io, "Identity to reference low at ", percent, " %")
end

"The segment is too short"
struct ErrorTooShort <: SegmentError
    len::UInt32
end

function Base.print(io::IO, x::ErrorTooShort)
    print(io, "Sequence too short at ", x,len, (isone(x.len) ? " base" : " bases"))
end

"Too many bases are insignificantly called in the sequence"
struct ErrorInsignificant <: SegmentError
    n_insignificant::UInt32
end

function Base.print(io::IO, x::ErrorInsignificant)
    print(io,
        "Sequence has ", x.n_insignificant, " insignificant ",
        (isone(x.n_insignificant) ? "base" : "bases")
    )
end

"Too many bases or amino acids are ambiguous"
struct ErrorAmbiguous <: SegmentError
    n_ambiguous::UInt32
end

function Base.print(io::IO, x::ErrorAmbiguous)
    print(io,
        "Sequence has ", x.n_ambiguous, " ambiguous ",
        (isone(x.n_ambiguous) ? "base" : "bases")
    )
end

"Some particular bases have too low depth"
struct ErrorLowDepthBases <: SegmentError
    n::UInt32
end

function Base.print(io::IO, x::ErrorLowDepthBases)
    print(io, "Sequence has ", x.n, " low-depth ", (isone(x.n) ? "base" : "bases"))
end

"Fraction of reference covered by reads/query is too low"
struct ErrorLowCoverage <: SegmentError
    coverage::Float32
end

function Base.print(io::IO, x::ErrorLowCoverage)
    n = @sprintf("%.3f", x.coverage)
    print(io, "Coverage is low at ", n)
end

"N'th round of assembly is too different from N-1'th round"
struct ErrorAssemblyNotConverged <: SegmentError
    identity::Float32
end

function Base.print(io::IO, x::ErrorAssemblyNotConverged)
    percent = round(x.identity * 100, digits=1)
    print(io, "Assembly not converged, at ", percent, " % identity")
end

"Sequence is flanked by invalid sequences - probably linkers or primers"
struct ErrorLinkerContamination <: SegmentError
    fiveprime:: Union{Nothing, UInt32}
    threeprime::Union{Nothing, UInt32}

    function ErrorLinkerContamination(fiveprime, threeprime)
        fp = convert(Union{Nothing, UInt32}, fiveprime)
        tp = convert(Union{Nothing, UInt32}, threeprime)
        if fp === tp === nothing
            throw(ArgumentError("Both fields cannot be `nothing`"))
        end
        new(fp, tp)
    end
end

function Base.print(io::IO, x::ErrorLinkerContamination)
    s = "Linker/primer contamination at ends, check "
    both = (x.fiveprime !== nothing) & (x.threeprime !== nothing)
    if x.fiveprime !== nothing
        s *= "first " * string(x.fiveprime) * (isone(x.fiveprime) ? " base" : " bases")
    end
    if x.threeprime !== nothing
        both && (s *= " and ")
        s *= "last " * string(x.threeprime) * (isone(x.threeprime) ? " base" : " bases")
    end
    print(io, s)
end

"The segment is missing a non-auxiliary protein"
struct ErrorMissingProtein <: SegmentError
    protein::Protein
end

function Base.print(io::IO, x::ErrorMissingProtein)
    print(io, "Missing non-auxiliary protein: \"", x.protein, '\"')
end

"Frameshift mutation"
struct ErrorFrameShift <: ProteinError
    indel::Indel
end

function Base.print(io::IO, x::ErrorFrameShift)
    print(io, "Frameshift: ", indel_message(x.indel))
end

"An indel is too big to be biologically plausible, or needs special attention"
struct ErrorIndelTooBig <: ProteinError
    indel::Indel
end

function Base.print(io::IO, x::ErrorIndelTooBig)
    print(io, "Indel too big: ", indel_message(x.indel))
end

"5' end of protein is deleted. This rarely happens naturally, and merits special attention"
struct ErrorFivePrimeDeletion <: ProteinError
    indel::Indel
end

function Base.print(io::IO, x::ErrorFivePrimeDeletion)
    print(io,
        "Deletion of ", length(x.indel),
        (isone(length(x.indel)) ? " base" : " bases"),
        " at 5' end"
    )
end

"Frameshift or substitution added a stop codon too early compared to reference"
struct ErrorEarlyStop <: ProteinError
    # We can't necessarily have expected pos, because the sequence may simply
    # stop before that part that aligns to the expected stop, so we can't
    # look at the alignment and see where the segment ought to stop.
    observed_pos::UInt32
    expected_naa::UInt32
    observed_naa::UInt32
end

function Base.print(io::IO, x::ErrorEarlyStop)
    print(
        io,
        "Protein stops early at segment pos ", x.observed_pos,
        " after ", x.observed_naa, " aa, reference is ",
        x.expected_naa, " aa"
    )
end

"Stop codon is mutated, protein stops later than expected"
struct ErrorLateStop <: ProteinError
    expected_pos::UInt32
    observed_pos::UInt32
    expected_naa::UInt32
    observed_naa::UInt32
end

function Base.print(io::IO, x::ErrorLateStop)
    print(
        io,
        "Protein stops late at segment pos ", x.observed_pos,
        " after ", x.observed_naa, " aa, reference stops at ",
        x.expected_stop, " after ", x.expected_naa, " aa"
    )
end

"ORF runs over edge of DNA sequence"
struct ErrorNoStop <: ProteinError end

function Base.print(io::IO, x::ErrorNoStop)
    print(io, "No stop codon")
end

"Length of coding sequence is not divisible by 3."
struct ErrorCDSNotDivisible <: ProteinError
    len::UInt32

    function ErrorCDSNotDivisible(x)
        len = convert(UInt32, x)
        if iszero(len % 3)
            throw(ArgumentError("Length must not be divisible by 3"))
        end
        new(len)
    end
end

function Base.print(io::IO, x::ErrorCDSNotDivisible)
    print(io, "CDS has length ", x.len, ", not divisible by 3")
end