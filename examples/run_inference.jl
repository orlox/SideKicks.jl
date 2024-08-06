using SideKicks
using BenchmarkTools
using CairoMakie
using Distributions
using HDF5

# USER sets this! - Ensure this matches a directory in examples/ which has the appropriate files and filenames
system_id = "vfts243" 

## 

# Obtain observations and priors from external files
obs_string = open("examples/$system_id/observations.jl", "r") do file
        read(file, String)
end
obs = eval(Meta.parse(obs_string))


priors_string = open("examples/$system_id/priors.jl", "r") do file
    read(file, String)
end
priors = eval(Meta.parse(priors_string))

mcmc_cauchy, props_cauchy = SideKicks.createGeneralMCMCModel( observations=obs, priors=priors, likelihood=:Cauchy)

##

use_general_model = true
if use_general_model
    which_model  = :general
else
    which_model  = :simplified
end

mcmcStruct = SideKicks.RunKickMCMC(
        which_model = which_model,
        observations = obs,
        priors = priors,
        nuts_warmup_count = 200,
        nuts_acceptance_rate = 0.8,
        nsamples = 200,
        nchains = 8)

##

SideKicks.SaveResults("examples/$system_id/results.hdf5", mcmcStruct, obs_string, priors_string)
