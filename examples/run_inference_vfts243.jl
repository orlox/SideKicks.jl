#=
# Inference example for VFTS 243

In this example we try to match the observed properties of VFTS 243. We start by loading
up the SideKicks package, as well the Distributions package which we use to define the priors.
=#

using SideKicks
using Distributions

##
#=
We set up our observations using the @Observations macro. This macro will return a tuple containing
an instance of `Observation`, as well as a string with the code needed to initialize the `Observation` instance.
The reason for this is that when we save results, we also want to save the precise observations that
were used for the MCMC.

Currently, The likelihoods are taken to be Normal distributions. This can potentially be generalized in
the future to arbitrary user-defined likelihood functions.
=#

obs = SideKicks.@Observations([
    [:P_f,  10.4031, 0.01,   day],
    [:e_f,  0.017,   0.012,  1],
    [:m1_f, 25.0,    2.3,    m_sun],
    [:K1,   81.4,    1.3,    km_per_s],
    #[:v_N,  143,     12,     km_per_s],
    #[:v_E,  408,     8,      km_per_s],
    [:v_N,  138.8,    7.6,     km_per_s], # Gaia 4' w/ MCMC
    [:v_E,  409.3,    9.3,     km_per_s], # Gaia 4' w/ MCMC
    [:v_r,  261.5,   0.42,    km_per_s], # Almeida
    [:ω_f,  66,      53,     degree]
]) 

##
#=
Next, we choose our priors. In here we are using distributions defined in the `Distributions.jl`` package.
=#
priors = SideKicks.@Priors(
    logm1_dist = Uniform(0.1,3), # in log(Msun)
    logm2_dist = Uniform(0.1,3), # in log(Msun)
    logP_dist  = Uniform(-1,3),   # in log(days)
    vkick_dist = Uniform(0,4), # in 100 km/s
    frac_dist  = Uniform(0,1.0),
    e_dist = Uniform(0,0.01),
    #venv_N_100kms_dist = Normal(145/100, 12/100),
    #venv_E_100kms_dist = Normal(396/100, 12/100),
    #venv_r_100kms_dist = Normal(271.6/100, 12.2/100)
    venv_N_100kms_dist = Normal(146/100, 40/100), # Gaia 4' w/ MCMC
    venv_E_100kms_dist = Normal(393/100, 42/100), # Gaia 4' w/ MCMC
    venv_r_100kms_dist = Normal(270.3/100, 11.1/100) # Almeida w/ MCMC
)

##
#=
Finally, we run the MCMC using the NUTS algorithm. In the example below we are computing 3 chains with 100 samples each.
This in practice is very low, and a ballpark suggestion would be to total at least $\sim10^4$ samples for any meaningful
analysis, and an order of magnitude higher for final science runs.
=#

kick_mcmc = SideKicks.KickMCMC(
        which_model = :general,
        observations = obs,
        priors = priors,
        nuts_warmup_count = 100,
        nuts_acceptance_rate = 0.8,
        nsamples = 100,
        nchains = 4)

##
#=
Results from the MCMC can be saved to an HDF5 file. In practice, it is ideal to separate the computation of the MCMC
and its analysis, so plotting of the results is done separetely.
=#
SideKicks.SaveResults("vfts243_results.hdf5", kick_mcmc)
