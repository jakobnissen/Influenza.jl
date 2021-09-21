"""
This module contains generally useful functionality related to influenza processing. It is imported by other packages.

It is distinct from `InfluenzaCore` in that it is a more heavyweight dependency, and distinct from other packages in that this package does not serve any particular purpose, but as a swiss knife tool for Influenza analysis.
"""
module Influenza

using ErrorTypes: Option, some, none, is_error, unwrap, unwrap_or, @unwrap_or, and_then
using FASTX: FASTA
using InfluenzaCore: Segment, Protein, Segments, SeroType, Proteins, source
using StructTypes: StructTypes # for JSON (de)serialization
using JSON3: JSON3
using BioSequences: BioSequences, LongDNASeq, LongAminoAcidSeq, NucleotideSeq,
    BioSequence, DNACodon, @mer_str, isgap, @aa_str, @biore_str, NucleicAcidAlphabet,
    DNA, DNA_Gap
using BioAlignments: BioAlignments
using Printf: Printf, @sprintf

const BA = BioAlignments

imap(f) = x -> Iterators.map(f, x)
ifilter(f) = x -> Iterators.filter(f, x)

"Conserved sequence at influenza segment 5' end"
const TERMINAL_INFLUENZA_5 = mer"AGCAAAAGCAGG"dna

"Conserved sequence at influenza segment 3' end"
const TERMINAL_INFLUENZA_3 = mer"CTTGTTTCTCCT"dna

const INFLUENZA_VERSION = let
    p = pathof(Influenza)::String
    project = joinpath(dirname(dirname(p)), "Project.toml")
    toml = read(project, String)
    m = match(r"(*ANYCRLF)^version\s*=\s\"(.*)\"$"m, toml)::RegexMatch
    VersionNumber(m[1]::SubString{String})
end

"Singleton struct used for three-valued logic, in e.g. `Union{Bool, Maybe}`"
struct Maybe end

include("errors.jl")
include("alignment.jl")
include("assembly.jl")
include("serialization.jl")
include("blast.jl")

# This is type piracy, but hopefully it's OK!
StructTypes.StructType(::Type{<:BioSequence}) = StructTypes.StringType()
StructTypes.StructType(::Type{<:AbstractUnitRange}) = StructTypes.Struct()

export TERMINAL_INFLUENZA_5,
    TERMINAL_INFLUENZA_3,

    Assembly,
    AssemblyProtein,
    Reference,
    AlignedAssembly,

    translate_reference,

    translate_proteins,
    alignment_identity,
    ha0_cleavage,
    store_references,
    load_references,
    annotate,

    # Exports from InfluenzaCore
    Segment, Segments, SeroType, Proteins, Protein, source

end # module
