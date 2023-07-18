module aOCS
using RCall, DataFrames, Distributions

include("optisel.jl")
include("merge-split.jl")
include("reproduce.jl")
include("xps.jl")
include("sum.jl")

end # module aOCS
