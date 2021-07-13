"""
This module contains generally useful functionality related to influenza processing. It is imported by other packages.

It is distinct from `InfluenzaCore` in that it is a more heavyweight dependency, and distinct from other packages in that this package does not serve any particular purpose, but as a swiss knife tool for Influenza analysis.
"""
module Influenza

using InfluenzaCore
using BioSequences
using BioAlignments
using ErrorTypes

"Singleton struct used for three-valued logic, in e.g. Union{Bool, Maybe}"
struct Maybe end

include("alignment.jl")

export is_stop,
    alignment_identity,
    ha0_cleavage,
    DEFAULT_DNA_ALN_MODEL,
    DEFAULT_AA_ALN_MODEL,

    # Exports from InfluenzaCore
    Segment, Segments, SeroType, Proteins, Protein, source

end # module
