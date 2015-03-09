# TODO: make this into a real package
# TODO: some docs
# TODO: Integrate with TravisCI & Coveralls

module nbody

using OctTrees
import OctTrees: modify, map, stop_cond, getx, gety, getz

using GeometricalPredicates

# TODO: change these to submodules...

include("cosmology.jl")
include("space/common.jl")
include("space/newtonian.jl")
include("space/cosmological.jl")
include("particle.jl")
include("world/world.jl")
include("world/realizations.jl")
include("parallel.jl")
include("simulation.jl")
include("forces/common.jl")
include("forces/newtonian.jl")
include("forces/cosmological.jl")
include("optimization/optimization.jl")
include("optimization/gradients.jl")
include("plotting.jl")

#include("../test/test.jl")

end # module
