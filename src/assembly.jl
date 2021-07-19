abstract type InfluenzaError end
abstract type SegmentError <: InfluenzaError end
abstract type ProteinError <: InfluenzaError end

"Protein or DNA sequence has too low identity compared to the reference"
struct ErrorLowIdentity <: InfluenzaError
    identity::Float32
end

"Too many bases are insignificantly called in the sequence"
struct ErrorInsignificant <: SegmentError
    n_insignificant::UInt32
end

"Too many bases or amino acids are ambiguous"
struct ErrorAmbiguous <: SegmentError
    n_ambiguous::UInt32
end

"Mean depth across sequence is too low"
struct ErrorLowDepth <: SegmentError
    depth::Float32
end

"Fraction of reference covered by reads/query is too low"
struct ErrorLowCoverage <: SegmentError
    coverage::Float32
end

"N'th round of assembly is too different from N-1'th round"
struct ErrorAssemblyNotConverged <: SegmentError
    identity::Float32
end

"Sequence is flanked by invalid sequences - probably linkers or primers"
struct ErrorLinkerContamination <: SegmentError
    fiveprime:: Union{Nothing, UnitRange{UInt32}}
    threeprime::Union{Nothing, UnitRange{UInt32}}
end

"The segment is missing a non-auxiliary protein"
struct ErrorMissingProtein <: SegmentError
    protein::Protein
end

# Protein errors
"Struct representing an indel compared to the reference"
struct Indel
    # range: Position in ref (del) or asm (ins) of bases affected
    range::UnitRange{UInt32}
    # position: seq aligns between pos and pos+1 in the other seq
    position::UInt32
    is_deletion::Bool
end

Base.length(x::Indel) = length(x.range)

struct ErrorFrameShift <: ProteinError
    indel::Indel
end

"Too many indels in sequence - alignment probably went wrong"
struct ErrorTooManyIndels <: ProteinError
    indels::Vector{Indel}
end

"An indel is too big to be biologically plausible, or needs special attention"
struct ErrorIndelTooBig <: ProteinError
    indel::Indel
end

"5' end of protein is deleted. This rarely happens naturally, and merits special attention"
struct ErrorFivePrimeDeletion <: ProteinError
    indel::Indel
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

"Stop codon is mutated, protein stops later than expected"
struct ErrorLateStop <: ProteinError
    expected_pos::UInt32
    observed_pos::UInt32
    expected_naa::UInt32
    observed_naa::UInt32
end

"ORF runs over edge of DNA sequence"
struct ErrorNoStop <: ProteinError end

"Length of coding sequence is not divisible by 3."
struct ErrorCDSNotDivisible <: ProteinError
    len::UInt32
end

"""
ReferenceProtein

Struct that holds data of one ORF in a segment.
"""
struct ReferenceProtein
    var::Protein
    orfs::Vector{UnitRange{UInt32}}
end

# This constructor validates orfs - may not be necessary
function ReferenceProtein(
    protein::Protein,
    orfs::Vector{<:UnitRange{<:Unsigned}},
    seq::NucleotideSeq
)
    issorted(orfs) || sort!(orfs)
    seqlen = length(seq)
    for orf in orfs
        if isempty(orf) || iszero(first(orf)) || last(orf) > seqlen
            throw(BoundsError(seq, orf))
        end
    end
    ReferenceProtein(protein, orfs)
end

"""
    A Reference that an assembly can be compared against
"""
struct Reference
    name::String
    segment::Segment
    seq::LongDNASeq
    proteins::Vector{ReferenceProtein}
end

struct Assembly
    name::String
    # none means unknown
    segment::Option{Segment}
    # none means no bases are insignificant, or unknown
    insignificant::Option{BitVector}
    seq::LongDNASeq
end

function Assembly(record::FASTA.Record, segment::Union{Segment, Nothing}, check_significance::Bool=true)
    itr = (i in UInt8('a'):UInt8('z') for i in @view record.data[record.sequence])
    insignificant = check_significance && any(itr) ? some(BitVector(itr)) : none(BitVector)
    name = let
        header = FASTA.header(rec)
        header === nothing ? "" : header
    end
    seq = FASTA.sequence(LongDNASeq, record)
    sgmt = segment === nothing ? none(Segment) : some(segment)
    return Assembly(name, sgmt, insignificant, seq)
end

# TODO: Should the errors be part of the struct, or stored outside the struct?
# they are not intrinsic to the information in it, but derived from it.
"""
    AssemblyProtein

A struct to store the information about a protein in an assembly, which has
been compared to its reference.
"""
struct AssemblyProtein
    variant::Protein
    orfs::Vector{UnitRange{UInt32}}
    identity::Float64
    errors::Vector{ProteinError}
end

function AssemblyProtein(
    protein::ReferenceProtein,
    aln::PairwiseAlignment{LongDNASeq, LongDNASeq},
    ref::Reference
)
    coding_mask = falses(length(ref.seq))
    for orf in protein.orfs
        coding_mask[orf] .= true
    end

    orfseq, orfs, errors, indels = compare_proteins_in_alignment(protein, coding_mask, aln)
    aaseq = BioSequences.translate(orfseq)
    ref_aas = (i for (i,n) in zip(ref.seq, coding_mask) if n)
    # Last 3 nts are the stop codon
    refaa = BioSequences.translate(LongDNASeq(collect(ref_aas)[1:end-3]))
    aaaln = pairalign(GlobalAlignment(), aaseq, refaa, DEFAULT_AA_ALN_MODEL).aln
    @assert aaaln !== nothing
    identity = alignment_identity(aaaln)::Float64
    return AssemblyProtein(protein.var, orfs, identity, errors)
end

is_stop(x::DNACodon) = (x === mer"TAA") | (x === mer"TAG") | (x === mer"TGA")

"Adds one nt at the end of the codon, moving it. If nt is ambiguous, return `nothing`"
function push_codon(x::DNACodon, nt::DNA)
    val = @inbounds BioSequences.twobitnucs[reinterpret(UInt8, nt) + 1]
    enc = (reinterpret(UInt64, x) << 2 | val) & UInt64(0x3f)
    ifelse(val === 0xff, nothing, reinterpret(DNACodon, enc))
end

"""Compares a protein and an alignment between a segment containing the protein
and the referece segment that contains that protein.
"""
function compare_proteins_in_alignment(
    protein::ReferenceProtein,
    coding_mask::BitVector,
    aln::PairwiseAlignment{LongDNASeq, LongDNASeq}
)::Tuple{LongDNASeq, Vector{UnitRange{UInt32}}, Vector{ProteinError}, Vector{Indel}}
    nucleotides = sizehint!(DNA[], 1200)
    last_coding_ref_pos = last(last(protein.orfs))
    codon = mer"AAA" # arbitrary starting codon
    seg_pos = ref_pos = n_deletions = n_insertions = 0
    fiveprime_truncated = 0
    maybe_expected_stop = none(Int)
    errors = ProteinError[]
    indels = Indel[]
    orfs = UnitRange{UInt32}[]
    seg_orfstart = nothing


    for (seg_nt, ref_nt) in aln
        seg_pos += (seg_nt !== DNA_Gap)
        ref_pos += (ref_nt !== DNA_Gap)
        is_coding = !iszero(ref_pos) && coding_mask[ref_pos]

        # Check for 5' truncation
        if iszero(seg_pos)
            fiveprime_truncated += is_coding
        else
            if !iszero(fiveprime_truncated)
                indel = Indel(
                    UInt32(ref_pos - n_deletions):UInt32(ref_pos - 1),
                    0,
                    true
                )
                push!(errors, ErrorFivePrimeDeletion(indel))
            end
            fiveprime_truncated = 0
        end

        # Add ORF if is coding and update seg_orfstart if applicable
        if is_coding
            if (seg_orfstart === nothing) & (seg_nt !== DNA_Gap)
                seg_orfstart = seg_pos
            end
        else
            if seg_orfstart !== nothing
                push!(orfs, UInt32(seg_orfstart):UInt32(seg_pos - 1))
                seg_orfstart = nothing
            end
        end

        if ref_pos == last_coding_ref_pos
            maybe_expected_stop = some(Int(seg_pos))
        end

        # All the rest of the operations only make sense if
        # the sequence is coding
        is_coding || continue

        # Check for deletions and update the codon
        if seg_nt == DNA_Gap
            n_deletions += 1
        else
            codon = let
                p = push_codon(codon, seg_nt)
                p === nothing ? mer"AAA" : p
            end
            push!(nucleotides, seg_nt)
            if !iszero(n_deletions)
                indel = Indel(
                    UInt32(ref_pos - n_deletions):UInt32(ref_pos - 1),
                    seg_pos - 1,
                    true
                )
                push!(indels, indel)
                if !iszero(length(indel) % 3)
                    push!(errors, ErrorFrameShift(indel))
                end
                if length(indel) > 21
                    push!(errors, ErrorIndelTooBig(indel))
                end
                n_deletions = 0
            end
        end

        # Check for insertions
        if ref_nt == DNA_Gap
            n_insertions += 1
        elseif !iszero(n_insertions)
            indel = Indel(
                UInt32(seg_pos - n_deletions):UInt32(seg_pos - 1),
                ref_pos - 1,
                false
            )
            push!(indels, indel)
            if !iszero(length(indel) % 3)
                push!(errors, ErrorFrameShift(indel))
            end
            if length(indel) > 36
                push!(errors, ErrorIndelTooBig(indel))
            end
            n_insertions = 0
        end

        # Only stop if we find a stop codon NOT in an intron
        if is_stop(codon) && iszero(length(nucleotides) % 3)
            n_aa = div(length(nucleotides), 3)
            expected_n_aa = div(sum(coding_mask), 3)

            # If we haven't yet reached the point where the stop ought to be
            if is_error(maybe_expected_stop)
                push!(errors, ErrorEarlyStop(seg_pos, expected_n_aa, n_aa))
            else
                expected_stop = unwrap(maybe_expected_stop)
                if expected_stop != seg_pos
                    @assert seg_pos > expected_stop
                    push!(errors, ErrorLateStop(expected_stop, seg_pos, expected_n_aa, n_aa))
                end
            end
            break
        end
    end

    # Add final orf after loop
    if seg_orfstart !== nothing
        push!(orfs, UInt16(seg_orfstart):UInt16(seg_pos))
    end

    dnaseq = LongDNASeq(nucleotides)
    # Is seq length divisible by three?
    remnant = length(dnaseq) % 3
    if !iszero(remnant)
        push!(errors, ErrorCDSNotDivisible(length(dnaseq)))
    end

    # Does it end with a stop? If so, remove it, else report error
    dnaseq = iszero(remnant) ? dnaseq : dnaseq[1:end-remnant]
    if isempty(dnaseq) || !is_stop(DNACodon(dnaseq[end-2:end]))
        push!(errors, ErrorNoStop())
    else
        dnaseq = dnaseq[1:end-3]
    end

    # If there are too many indel messages, that's presumably a problem
    if length(indels) > 3
        push!(errors, ErrorTooManyIndels(indels))
    end

    return dnaseq, orfs, errors, indels
end


"""
    AlignedAssembly

Struct to store information about a DNA sequence aligned to its reference.
"""
struct AlignedAssembly
    assembly::Assembly
    reference::Reference
    aln::PairwiseAlignment{LongDNASeq, LongDNASeq}
    identity::Float64
    proteins::Vector{AssemblyProtein}
    errors::Vector{Union{ErrorLowIdentity, SegmentError}}
end

function AlignedAssembly(asm::Assembly, ref::Reference)
    if unwrap_or(asm.segment, ref.segment) !== ref.segment
        error("Cannot make AlignedAssembly of different segments")
    end

    # For optimization: The large majority of time is spent on this alignment
    aln = pairalign(OverlapAlignment(), asm.seq, ref.seq, DEFAULT_DNA_ALN_MODEL).aln
    @assert aln !== nothing

    identity = alignment_identity(aln)::Float64

    proteins = map(ref.proteins) do protein
        AssemblyProtein(protein, aln, ref)
    end

    errors = Union{ErrorLowIdentity, SegmentError}[]
    
    # Low identity to reference
    identity < 0.9 && push!(errors, ErrorLowIdentity(identity))

    # Insignificant bases
    n_insignificant = unwrap_or(and_then(count, Int, asm.insignificant), 0)
    iszero(n_insignificant) || push!(errors, ErrorInsignificant(n_insignificant))

    # Ambiguous bases
    n_amb = count(isambiguous, asm.seq)
    iszero(n_amb) || push!(errors, ErrorAmbiguous(n_amb))

    return AlignedAssembly(asm, ref, aln, identity, proteins, errors)
end

