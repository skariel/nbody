module nbody

using OctTrees
import OctTrees: modify, map, stop_cond, getx, gety, getz

using GeometricalPredicates

include("types.jl")
include("simulation.jl")
include("compiletree.jl")
include("forces.jl")
include("realizations.jl")
include("plotting.jl")
include("backdynamics.jl")
include("test.jl")

end # module
