module Influenza

using InfluenzaCore
using BioSequences
using BioAlignments
using ErrorTypes

include("alignment.jl")
include("interactive.jl")

export is_stop,
    alignment_identity,
    DEFAULT_DNA_ALN_MODEL,
    DEFAULT_AA_ALN_MODEL,

    # Exports from InfluenzaCore
    Segment, Segments, SubType, SubTypes, Proteins, Protein, source

end # module
