"""
    DEFAULT_DNA_ALN_MODEL

Affine Gap Score model with parameters empirically chosen to strike a sensible balance between
indels and substitutions for DNA alignments.
"""
const DEFAULT_DNA_ALN_MODEL = AffineGapScoreModel(EDNAFULL, gap_open=-25, gap_extend=-2)

"""
    DEFAULT_AA_ALN_MODEL

Affine Gap Score model with parameters empirically chosen to strike a sensible balance between
indels and substitutions for amino acid alignments.
"""
const DEFAULT_AA_ALN_MODEL = AffineGapScoreModel(BLOSUM62, gap_open=-10, gap_extend=-2)

"""
    alignment_identity([m::AbstractAlignment], ::PairwiseAlignment)

Calculate the alignment identity between two sequences. Alignment is calculated as
n_matches divided by the ungapped length of the longest seq.
If `m` is `GlobalAlignment()`, as is default, include all nucleotides in alignment.
If `m` is `OverlapAlignment()`, do not include positions with gaps at the end.

If the shorter seq has zero length, returns `nothing`.
"""
function alignment_identity end

function alignment_identity(aln::PairwiseAlignment{T, T}) where {T <: BioSequence}
    alignment_identity(GlobalAlignment(), aln)
end

function alignment_identity(::GlobalAlignment, aln::PairwiseAlignment{T, T}) where {T <: BioSequence}
    alignment_identity(collect(aln))
end

function alignment_identity(::OverlapAlignment, aln::PairwiseAlignment{T, T}) where {T <: BioSequence}
    v = collect(aln)
    p = findfirst(t -> !isgap(t[1]) & !isgap(t[2]), v)
    p === nothing && return nothing
    pend = findlast(t -> !isgap(t[1]) & !isgap(t[2]), v)::Int
    alignment_identity(view(v, p:pend))
end

function alignment_identity(v::AbstractVector{Tuple{S, S}}) where {S <: BioSequences.BioSymbol}
    n_ident = len_query = len_subject = 0
    for (seqnt, refnt) in v
        n_ident += seqnt == refnt
        len_query += !isgap(seqnt)
        len_subject += !isgap(refnt)
    end
    len_longest = max(len_query, len_subject)
    iszero(len_longest) && return nothing
    return n_ident / len_longest
end


"""
    ha0_cleavage(::LongAminoAcidSeq)

Detects and returns the HA0 cleavage site. Returns a tuple `(site, is_hpai)`, where:
* `site` is a `LongAminoAcidSeq` if a site was found, otherwise `nothing`.
* `is_hpai` is `nothing` if `site` is, `Maybe()` if the pathogenicity cannot be determined,
`true` if it is HPAI, and `false` if it is LPAI.

# Examples
```julia-repl
julia> ha0_cleavage(aa"LATGLRNSPLREKRRKRGLFGAIAGFIEGGW")
(PLREKRRKRGLF, true)
```
"""
function ha0_cleavage(
    seq::LongAminoAcidSeq
)::Tuple{Union{Nothing, LongAminoAcidSeq}, Union{Nothing, Maybe, Bool}}
    motif = cleavage_motif(seq)
    # A nothing here means the motif was not properly detected
    motif === nothing && return (nothing, nothing)
    seq, pos = motif 
    return seq, is_leading_hpai(seq[pos-5:pos-1])
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
