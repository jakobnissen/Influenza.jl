const DEFAULT_DNA_ALN_MODEL = AffineGapScoreModel(EDNAFULL, gap_open=-25, gap_extend=-2)
const DEFAULT_AA_ALN_MODEL = AffineGapScoreModel(BLOSUM62, gap_open=-10, gap_extend=-2)

is_stop(x::DNACodon) = (x === mer"TAA") | (x === mer"TAG") | (x === mer"TGA")

"Alignment is calulated as n_matches divided by the length of the shortest seq"
function alignment_identity(aln::PairwiseAlignment{T, T})::Option{Float64} where {T <: BioSequence}
    n_ident = len_query = len_subject = 0
    for (seqnt, refnt) in aln
        n_ident += seqnt == refnt
        len_query += !is_gap(seqnt)
        len_subject += !is_gap(refnt)
    end
    len_smallest = min(len_query, len_subject)
    iszero(len_smallest) && return none
    return some(n_ident / len_smallest)
end