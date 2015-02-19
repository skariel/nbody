module nbody

using OctTrees
import OctTrees: modify, stop_cond, getx, gety, getz

using GeometricalPredicates

include("types.jl")
include("simulation.jl")
include("forces.jl")
include("realizations.jl")
include("plotting.jl")
include("backdynamics.jl")
include("compiletree.jl")
include("test.jl")

end # module
