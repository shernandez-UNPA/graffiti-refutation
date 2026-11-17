module GraffitiDE

include("invariants.jl")
include("conjectures.jl")
include("decoders.jl")
include("search.jl")

export CONJECTURES, Conjecture, violation, Params, Result, search, verify,
       decode, Decoder, graph6, edgelist, dist_matrix, spec, harmonic, randic,
       gravity_matrix, TOL

end
