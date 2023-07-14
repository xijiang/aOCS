module aOCS
using RCall, DataFrames, Distributions

include("optisel.jl")
include("merge-split.jl")
include("reproduce.jl")
include("xps.jl")

end # module aOCS
