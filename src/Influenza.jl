"""
This module contains generally useful functionality related to influenza processing. It is imported by other packages.

It is distinct from `InfluenzaCore` in that it is a more heavyweight dependency, and distinct from other packages in that this package does not serve any particular purpose, but as a swiss knife tool for Influenza analysis.
"""
module Influenza

using ErrorTypes
using FASTX
using InfluenzaCore
using BioSequences
using BioAlignments
using Serialization

const INFLUENZA_VERSION = let
    project = joinpath(dirname(dirname(pathof(Influenza))), "Project.toml")
    toml = read(project, String)
    m = match(r"(*ANYCRLF)^version\s*=\s\"(.*)\"$"m, toml)
    VersionNumber(m[1])
end

"Singleton struct used for three-valued logic, in e.g. Union{Bool, Maybe}"
struct Maybe end

include("assembly.jl")
include("alignment.jl")
include("serialization.jl")

export is_stop,
    alignment_identity,
    ha0_cleavage,
    DEFAULT_DNA_ALN_MODEL,
    DEFAULT_AA_ALN_MODEL,
    store_references,
    load_references,

    # Exports from InfluenzaCore
    Segment, Segments, SeroType, Proteins, Protein, source

end # module
