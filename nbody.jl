module nbody

using OctTrees
import OctTrees: modify, stop_cond, getx, gety, getz

include("types.jl")
include("simulation.jl")
include("forces.jl")
include("realizations.jl")
include("test.jl")
include("plotting.jl")

end # module
