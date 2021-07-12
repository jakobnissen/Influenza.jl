module Influenza

using InfluenzaCore
using BioSequences
using BioAlignments
using ErrorTypes

"Used for three-valued logic, in e.g. Union{Bool, Maybe}"
struct Maybe end

include("alignment.jl")

export is_stop,
    alignment_identity,
    cleavage,
    DEFAULT_DNA_ALN_MODEL,
    DEFAULT_AA_ALN_MODEL,

    # Exports from InfluenzaCore
    Segment, Segments, SeroType, Proteins, Protein, source

end # module
