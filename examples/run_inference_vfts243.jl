using SideKicks
using BenchmarkTools
using Distributions

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

kick_mcmc = SideKicks.KickMCMC(
        which_model = :general,
        observations = obs,
        priors = priors,
        nuts_warmup_count = 200,
        nuts_acceptance_rate = 0.8,
        nsamples = 200,
        nchains = 8)

##

SideKicks.SaveResults(String(@__DIR__)*"/vfts243/results.hdf5", kick_mcmc)
