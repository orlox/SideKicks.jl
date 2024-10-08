"""
Main module for `SideKicks.jl` -- a analysis package for performing parameter inference on compact object stellar binaries.

"""
module SideKicks

include("Constants.jl")
include("Orbits.jl")
include("BHFallback.jl")
include("Observations.jl")
include("Priors.jl")
include("TuringModels.jl")
include("KickMCMC.jl")
include("CornerPlots.jl")
include("OutputStorage.jl")

end # module
