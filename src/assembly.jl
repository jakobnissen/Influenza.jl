"""
ReferenceProtein

Struct that holds data of one ORF in a segment.
"""
struct ReferenceProtein
    var::Protein
    orfs::Vector{UnitRange{UInt32}}
end
StructTypes.StructType(::Type{ReferenceProtein}) = StructTypes.Struct()

# This constructor validates orfs - may not be necessary
function ReferenceProtein(
    protein::Protein,
    orfs::Vector{<:UnitRange{<:Unsigned}},
    seq::NucleotideSeq,
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
A Reference that an assembly can be compared against.

A reference holds a name, a segment, a DNA sequence, and a vector of `ReferenceProtein`, which gives the proteins encoded by the segment and their open reading frames.
"""
struct Reference
    name::String
    segment::Segment
    seq::LongDNASeq
    proteins::Vector{ReferenceProtein}
end
StructTypes.StructType(::Type{Reference}) = StructTypes.Struct()

"""
A DNA sequence representing an influenza segment.

Assemblies consists of a name, a DNA sequence, and optionally a `Segment` and a
bitvector, signifying the bases that are insignificantly called.

# Examples
```
julia> asm = Assembly("myseq", dna"ACC")
Assembly("myseq", ACC, none(Segment), none(BitVector))

julia> asm = Assembly("myseq2", dna"TC", some(Segments.PB1), some(trues(3)))
Assembly("myseq2", TC, some(InfluenzaCore.Segments.PB1), some(Bool[1, 1, 1]))
```
"""
struct Assembly
    name::String
    seq::LongDNASeq
    # none means unknown
    segment::Option{Segment}
    # none means no bases are insignificant, or unknown
    insignificant::Option{BitVector}
end

function Assembly(name::AbstractString, seq::BioSequence{<:NucleicAcidAlphabet})
    return Assembly(
        convert(String, name),
        convert(LongDNASeq, seq),
        none(Segment),
        none(BitVector),
    )
end

function Assembly(
    record::FASTA.Record,
    segment::Union{Segment, Nothing},
    check_significance::Bool=true,
)
    itr = (i in UInt8('a'):UInt8('z') for i in @view record.data[record.sequence])
    insignificant = check_significance && any(itr) ? some(BitVector(itr)) : none(BitVector)
    name = let
        header = FASTA.header(record)
        header === nothing ? "" : header
    end
    seq = FASTA.sequence(LongDNASeq, record)
    sgmt = segment === nothing ? none(Segment) : some(segment)
    return Assembly(name, seq, sgmt, insignificant)
end

"""
A struct to store the information about a protein in an assembly, which has
been compared to its reference. See the fields of the struct for its information.
"""
struct AssemblyProtein
    variant::Protein
    orfs::Option{Vector{UnitRange{UInt32}}}
    identity::Option{Float64}
    errors::Vector{ProteinError}
end

function AssemblyProtein(
    protein::ReferenceProtein,
    aln::BA.PairwiseAlignment{LongDNASeq, LongDNASeq},
    ref::Reference,
)
    coding_mask = falses(length(ref.seq))
    for orf in protein.orfs
        coding_mask[orf] .= true
    end

    orfseq, orfs, errors, indels = compare_proteins_in_alignment(protein, coding_mask, aln)
    remainder = rem(length(orfseq), 3)
    iszero(remainder) || resize!(orfseq, length(orfseq) - remainder)
    aaseq = BioSequences.translate(orfseq)
    ref_aas = (i for (i, n) in zip(ref.seq, coding_mask) if n)
    # Last 3 nts are the stop codon
    refaa = BioSequences.translate(LongDNASeq(collect(ref_aas)[1:(end - 3)]))
    aaaln = BA.pairalign(BA.GlobalAlignment(), aaseq, refaa, DEFAULT_AA_ALN_MODEL).aln
    @assert aaaln !== nothing

    # The orfseq can be empty if the alignment has sufficiently low identity.
    # in this case, we will store orfs and identity as none.
    (identity, orfs) = if isempty(orfseq)
        none(Float64), none(Vector{UnitRange{UInt32}})
    else
        some(alignment_identity(BA.GlobalAlignment(), aaaln)::Float64), some(orfs)
    end
    return AssemblyProtein(protein.var, orfs, identity, errors)
end

# I absolutely hate doing this. But if there is a deletion exactly around the edges of
# an orf - which can happen spuriously if the exact position of an indel is unknowable,
# then the AlignedAssembly will not work because the ORF positions will be thrown off.
# This function moves a deletion in the query a few bases if so.
# E.g in the start codon ATG here:
#  seq:    1 GGAAAA--TGATTGCA
#            ||||||  ||||||||
#  ref:    1 GGAAAAAATGATTGCA
# It's clear the indel is misplaced.
function reorder_deletions!(aln::BA.PairwiseAlignment, p::Integer, start::Bool)
    anchors = aln.a.aln.anchors
    anchorpos =
        searchsortedfirst(anchors, p; by=a -> a isa BA.AlignmentAnchor ? a.refpos : a)
    anchor = anchors[anchorpos]
    # Whether we should move bases from before deletion or after deletion.
    # we move bases from outside the ORF.
    delta = ifelse(start, -1, 1)
    in(anchorpos + 2delta, eachindex(anchors)) || return nothing
    anchor.op == BA.OP_DELETE || return nothing
    n_move = anchor.refpos - p + 1
    # Let's not move too many bases
    n_move > 10 && return nothing
    next = anchors[anchorpos + delta]
    next.op ∈ (BA.OP_MATCH, BA.OP_SEQ_MISMATCH, BA.OP_SEQ_MATCH) || return nothing
    next2 = anchors[anchorpos + 2delta]
    # If there are not enough bases to move, we just return nothing.
    # it is possible to do it even if there are not enough bases, but it requires
    # this function to rewrite the number of anchors, so it complicates it.
    # This probably happens rarely if the sequences are otherwise good.
    if abs(next2.seqpos - next.seqpos) ≤ n_move
        return nothing
    end

    # Check if the modified alignment has as many matches as the old one
    query = if start
        aln.a.seq[(next.seqpos - n_move + 1):(next.seqpos)]
    else
        aln.a.seq[(anchor.seqpos + 1):(anchor.seqpos + n_move)]
    end
    ref1 = if start
        aln.b[(next.refpos - n_move + 1):(next.refpos)]
    else
        aln.b[(anchor.refpos + 1):(anchor.refpos + n_move)]
    end
    ref2 = aln.b[(p - n_move + 1):p]
    if count(Base.splat(==), zip(query, ref1)) > count(Base.splat(==), zip(query, ref2))
        return nothing
    end

    anchors[anchorpos] = BA.AlignmentAnchor(
        anchor.seqpos + n_move * delta,
        anchor.refpos + n_move * delta,
        BA.OP_DELETE,
    )
    anchors[anchorpos + delta] = BA.AlignmentAnchor(
        next.seqpos + n_move * delta,
        next.refpos + n_move * delta,
        next.op,
    )
    return nothing
end

function reorder_deletions!(aln::BA.PairwiseAlignment, orf::AbstractUnitRange)
    reorder_deletions!(aln, first(orf), true)
    reorder_deletions!(aln, last(orf), false)
end

"""
    is_stop(x::DNACodon)

Return whether the DNA Codon (a 3-mer) is TAA, TAG or TGA.
"""
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
    aln::BA.PairwiseAlignment{LongDNASeq, LongDNASeq},
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
        is_intron = !is_coding && ref_pos < last_coding_ref_pos

        # Check for 5' truncation
        if iszero(seg_pos)
            fiveprime_truncated += is_coding
        else
            if !iszero(fiveprime_truncated)
                indel = Indel(UInt32(ref_pos - n_deletions):UInt32(ref_pos - 1), 0, true)
                push!(errors, ErrorFivePrimeDeletion(indel))
            end
            fiveprime_truncated = 0
        end

        # Add ORF if not in intron and update seg_orfstart if applicable
        if is_intron
            if seg_orfstart !== nothing
                push!(orfs, UInt32(seg_orfstart):UInt32(seg_pos - 1))
                seg_orfstart = nothing
            end
        else
            if (seg_orfstart === nothing) & (seg_nt !== DNA_Gap)
                seg_orfstart = seg_pos
            end
        end

        # We expect to stop 3 bases after last coding position (at a STOP codon)
        if ref_pos == last_coding_ref_pos + 3
            maybe_expected_stop = some(Int(seg_pos))
        end

        # All the rest of the operations only make sense if
        # the sequence is coding. We skip the introns since they are spliced
        # out of the mRNA, but any pos after last intron must be terminated
        # by encountering a stop codon
        is_intron && continue

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
                    true,
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
                UInt32(seg_pos - n_insertions):UInt32(seg_pos - 1),
                ref_pos - 1,
                false,
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
                    push!(
                        errors,
                        ErrorLateStop(expected_stop, seg_pos, expected_n_aa, n_aa),
                    )
                end
            end
            break
        end
    end

    # If we ended due to a stop, remove the stop codon
    if is_stop(codon) && iszero(length(nucleotides) % 3)
        resize!(nucleotides, length(nucleotides) - 3)
        seg_pos -= 3
    else
        push!(errors, ErrorNoStop())
    end

    # Add final orf after loop.
    seg_orfstart === nothing || push!(orfs, UInt16(seg_orfstart):UInt16(seg_pos))
    dnaseq = LongDNASeq(nucleotides)
    return dnaseq, orfs, errors, indels
end

"""
    AlignedAssembly(::Assembly, ::Reference, force_termini=false)

Struct to store information about a DNA sequence aligned to its reference.
Creating this object automatically aligns the assembly to the reference and validates
it, adding any errors to its `errors` field.

If `force_termini` is set to true, then the conserved 5' and 3' motifs must be
present in the references, or it throws an error. Setting this enables any errors
related to these termini in the `Assembly` to be in the `errors` fields.

See the fields of the struct for the information contained.
"""
struct AlignedAssembly
    assembly::Assembly
    reference::Reference
    aln::BA.PairwiseAlignment{LongDNASeq, LongDNASeq}
    identity::Float64
    proteins::Vector{AssemblyProtein}
    errors::Vector{Union{ErrorLowIdentity, SegmentError}}
end

function AlignedAssembly(asm::Assembly, ref::Reference, force_termini::Bool=false)
    if unwrap_or(asm.segment, ref.segment) !== ref.segment
        error("Cannot make AlignedAssembly of different segments")
    end

    # For optimization: The large majority of time is spent on this alignment
    aln = BA.pairalign(BA.OverlapAlignment(), asm.seq, ref.seq, DEFAULT_DNA_ALN_MODEL).aln
    @assert aln !== nothing

    # Move around deletions in query if they straddle an edge of an ORF.
    # Ugly hack, but I don't know how to deal with it otherwise
    for protein in ref.proteins, orf in protein.orfs
        reorder_deletions!(aln, orf)
    end

    identity = let
        id = alignment_identity(BA.OverlapAlignment(), aln)
        id === nothing ? 0.0 : id
    end::Float64

    proteins = map(ref.proteins) do protein
        AssemblyProtein(protein, aln, ref)
    end

    errors = Union{ErrorLowIdentity, SegmentError}[]

    # Insignificant bases
    n_insignificant = unwrap_or(and_then(count, Int, asm.insignificant), 0)
    iszero(n_insignificant) || push!(errors, ErrorInsignificant(n_insignificant))

    # Ambiguous bases
    n_amb = count(BioSequences.isambiguous, asm.seq)
    iszero(n_amb) || push!(errors, ErrorAmbiguous(n_amb))

    # Termini. If we force termini to be present in ref, we add an error if it's
    # not present in the assembly.
    (m_notermini, m_linker) = check_termini(aln, force_termini, 2)
    if force_termini && m_notermini !== nothing
        push!(errors, m_notermini)
    end
    m_linker !== nothing && push!(errors, m_linker)

    return AlignedAssembly(asm, ref, aln, identity, proteins, errors)
end

function Reference(x::AlignedAssembly)
    proteins = map(x.proteins) do protein
        ReferenceProtein(protein.variant, unwrap(protein.orfs))
    end
    Reference(x.assembly.name, x.reference.segment, x.assembly.seq, proteins)
end

"""
    translate_reference(::AlignedAssembly)

Create a new `Reference` object from the assembly in `AlignedAssembly`. This function
does NOT perform any checks that the aligned assembly is suitable to make into a
reference, but may fail if the aligned assembly is so bad a sensible reference cannot
be constructed.
"""
function translate_reference(x::AlignedAssembly)
    proteins = map(x.proteins) do protein
        ReferenceProtein(protein.var, unwrap(protein.orfs))
    end
    return Reference(x.assembly.name, x.reference.segment, x.assembly.seq, proteins)
end

"""
    translate_proteins(::AlignedAssembly)

Get a vector of `Option{LongAminoAcidSeq}`, one from each protein of the aligned
assembly. If the length of the ORF is not divisible by 3, truncates bases from the 3' end.
Does not do any validation of the AA sequences.
"""
function translate_proteins(alnasm::AlignedAssembly)
    dnaseq = LongDNASeq()
    result = Option{LongAminoAcidSeq}[]
    for protein in alnasm.proteins
        if is_error(protein.orfs)
            push!(result, none(LongAminoAcidSeq))
        else
            empty!(dnaseq)
            for orf in unwrap(protein.orfs)
                append!(dnaseq, alnasm.assembly.seq[orf])
            end
            resize!(dnaseq, length(dnaseq) - length(dnaseq) % 3)
            push!(result, some(BioSequences.translate(dnaseq)))
        end
    end
    return result
end

"""
    check_termini(aln, force=true, max_subs=1) -> (ismissing, haslinker)

Check the alignment between a query and reference sequence starts and ends at the
correct positions. `ismissing` isa `Union{Nothing, ErrorNoTermini}`, and `haslinker`
isa `Union{Nothing, ErrorLinkerContamination}`.

The reference sequence must have the expected termini - this is checked by searching
for the known termini at the ends of the reference. Up to `max_subs` is allowed in this
check. If the termini are not found in the reference, `force` causes it to throw an error.
If `force` is `false`, return  `(nothing, nothing)`.
"""
function check_termini(
    aln::BA.PairwiseAlignment{<:NucleotideSeq, <:NucleotideSeq},
    force::Bool=true,
    max_subs::Integer=2,
)::Tuple{Union{Nothing, ErrorNoTermini}, Union{Nothing, ErrorLinkerContamination}}
    ref = aln.b
    rng_fw = BioSequences.approxsearch(
        ref,
        TERMINAL_INFLUENZA_5,
        max_subs,
        1,
        lastindex(TERMINAL_INFLUENZA_5) + max_subs,
    )
    rng_rv = BioSequences.approxrsearch(
        ref,
        TERMINAL_INFLUENZA_3,
        max_subs,
        lastindex(ref),
        lastindex(ref) - lastindex(TERMINAL_INFLUENZA_3) + 1 - max_subs,
    )
    if (isempty(rng_fw) | isempty(rng_rv))
        if force
            throw(ArgumentError("Terminal not found in reference"))
        else
            return (nothing, nothing)
        end
    end
    anchors = aln.a.aln.anchors
    isempty(anchors) && return (nothing, nothing)
    missfive, missthree = false, false
    contfive, contthree = nothing, nothing
    if anchors[2].op === BioAlignments.OP_DELETE
        missfive = true
    elseif anchors[2].op === BioAlignments.OP_INSERT
        contfive = UInt32(anchors[2].seqpos)
    end
    if anchors[end].op === BioAlignments.OP_DELETE
        missthree = true
    elseif anchors[end].op === BioAlignments.OP_INSERT
        contthree = UInt32(anchors[end].seqpos - anchors[end - 1].seqpos)
    end
    fst = !(missfive | missthree) ? nothing : ErrorNoTermini(missfive, missthree)
    scn = if (contfive === nothing) & (contthree === nothing)
        nothing
    else
        ErrorLinkerContamination(contfive, contthree)
    end
    return (fst, scn)
end
