"""
Main module for `SideKicks.jl` -- a analysis package for performing parameter inference on compact object stellar binaries.

"""
module SideKicks

include("Constants.jl")
include("Orbits.jl")
include("BHFallback.jl")
include("KickMCMC.jl")
include("KickViz.jl")

end # module
