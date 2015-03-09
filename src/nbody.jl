# TODO: make this into a real package
# TODO: some docs
# TODO: separate testing into its own directory
# TODO: Integrate with TravisCI & Coveralls

module nbody

using OctTrees
import OctTrees: modify, map, stop_cond, getx, gety, getz

using GeometricalPredicates

include("cosmology.jl")
include("types.jl")
include("parallel.jl")
include("simulation.jl")
include("forces/common.jl")
include("forces/newtonian.jl")
include("forces/cosmological.jl")
include("realizations.jl")
include("plotting.jl")
include("gradients.jl")
include("optimization.jl")

end # module
