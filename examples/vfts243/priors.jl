
priors = SideKicks.Priors(
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
