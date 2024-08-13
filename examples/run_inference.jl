using SideKicks
using BenchmarkTools
using Distributions

obs = SideKicks.@Observations([
    [:P,   10.4031,  0.01,   day],
    [:e,   0.017,    0.012,  1],
    [:m1,  25.0,     2.3,    m_sun],
    [:K1,  81.4,     1.3,    km_per_s],
    [:v_N,  143,     12,    km_per_s],
    [:v_E,  408,     8,     km_per_s],
    [:v_r,  260.2,   0.9,    km_per_s],
    [:ω,  66,   53,    degree]
]) 

priors = SideKicks.@Priors(
    logm1_dist = Uniform(0.1,3), # in log(Msun)
    logm2_dist = Uniform(0.1,3), # in log(Msun)
    logP_dist  = Uniform(-1,3),   # in log(days)
    vkick_dist = Exponential(1), # in 100 km/s
    frac_dist  = Uniform(0,1.0),
    e_dist = Uniform(0,0.01),
    v_N_100kms_dist = Normal(145/100, 12/100),
    v_E_100kms_dist = Normal(396/100, 12/100),
    v_r_100kms_dist = Normal(271.6/100, 12.2/100)
)

##

kick_mcmc = SideKicks.KickMCMC(
        which_model = :general,
        observations = obs,
        priors = priors,
        nuts_warmup_count = 200,
        nuts_acceptance_rate = 0.8,
        nsamples = 200,
        nchains = 8)

##

SideKicks.SaveResults("examples/vfts243_results.hdf5", kick_mcmc)