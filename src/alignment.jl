const DEFAULT_DNA_ALN_MODEL = AffineGapScoreModel(EDNAFULL, gap_open=-25, gap_extend=-2)
const DEFAULT_AA_ALN_MODEL = AffineGapScoreModel(BLOSUM62, gap_open=-10, gap_extend=-2)

is_stop(x::DNACodon) = (x === mer"TAA") | (x === mer"TAG") | (x === mer"TGA")

"Alignment is calulated as n_matches divided by the length of the shortest seq"
function alignment_identity(aln::PairwiseAlignment{T, T})::Option{Float64} where {T <: BioSequence}
    n_ident = len_query = len_subject = 0
    for (seqnt, refnt) in aln
        n_ident += seqnt == refnt
        len_query += !isgap(seqnt)
        len_subject += !isgap(refnt)
    end
    len_smallest = min(len_query, len_subject)
    iszero(len_smallest) && return none
    return some(n_ident / len_smallest)
end

"""
Check whether an amino acid sequence of HA is HPAI. Returns `nothing` if the cleavage
site was not detected, `Maybe()` if the HPAI-ness cannot easily be determined, and
a `Bool` otherwise.
"""
function is_hpai(seq::LongAminoAcidSeq)::Union{Nothing, Maybe, Bool}
    motif = cleavage_motif(seq)
    # A nothing here means the motif was not properly detected
    motif === nothing && return nothing
    seq, pos = motif 
    return is_leading_hpai(seq[pos-5:pos-1])
end

# Nothing if there is no motif, else (seq, pos_of_cleaving_R)
function cleavage_motif(seq::LongAminoAcidSeq)::Union{Nothing, Tuple{LongAminoAcidSeq, Int}}
    model = AffineGapScoreModel(BLOSUM62, gap_open=-10, gap_extend=-2)
    aln = pairalign(SemiGlobalAlignment(), seq, aa"LATGLRNSPLREKRRKRGLFGAIAGFIEGGW", model).aln
    site = AminoAcid[]
    refpos = 0
    for (q, r) in aln
        isgap(r) || (refpos += 1)
        refpos in 9:20 && !isgap(q) && push!(site, q)
    end
    motif = LongAminoAcidSeq(site)
    pos = findfirst(biore"[RK]GLF"aa, motif)
    
    # If the motif is not there, return nothing
    pos === nothing && return nothing

    # If the motif is not the last 4 AAs, return nothing
    last(pos) == length(motif) || return nothing

    # If the total site is too short, return nothing
    length(site) < 9 && return nothing
    return motif, first(pos)
end

# Returns nothing if it's uncertain
function is_leading_hpai(leading::LongAminoAcidSeq)::Union{Bool, Maybe}
    # Checks 5 bases before cleaving R/K residue. If 4/5 are basic, it's HPAI, if
    # <3, LPAI, if 3 it's ambiguous
    is_rk(aa) = (aa === AA_R) | (aa === AA_K)
    n_basic = sum(is_rk, leading)
    return if n_basic < 3
        false
    elseif n_basic == 3
        Maybe()
    else
        true
    end
end