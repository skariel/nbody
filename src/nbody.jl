# TODO: make this into a real package
# TODO: some docs
# TODO: Integrate with TravisCI & Coveralls

module nbody

using Compat

using OctTrees
import OctTrees: modify, map, stop_cond, getx, gety, getz

using GeometricalPredicates

# TODO: change these to submodules...

include("units.jl")
include("cosmology.jl")
include("space.jl")
include("particle.jl")
include("world/world.jl")
include("world/realizations.jl")
include("world/utils.jl")
include("parallel.jl")
include("simulation/common.jl")
include("simulation/newtonian.jl")
include("simulation/cosmological.jl")
include("forces/common.jl")
include("forces/newtonian.jl")
include("forces/cosmological.jl")
include("optimization/optimization.jl")
include("optimization/gradients.jl")
include("plotting.jl")
include("G2/IO.jl")
include("simulation/IO.jl")
include("world/IO.jl")

end # module
