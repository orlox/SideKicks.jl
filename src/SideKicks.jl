"""
Main module for `SideKicks.jl` -- a analysis package for performing parameter inference on compact object stellar binaries.

X functions are exported from this module for public use:

(RTW - how to find these?)
(RTW - how to replace Sidekicks.Sidekicks with just a single Sidekicks. This is done on the Documenter documentation page...)

"""
module SideKicks

include("Constants.jl")
include("Orbits.jl")
include("BHModels.jl")
include("KickMCMC.jl")
include("KickViz.jl")

end # module
