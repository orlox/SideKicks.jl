using SideKicks
using BenchmarkTools
using CairoMakie
using Distributions

obs = SideKicks.createObservations([
    [:P,   10.4031,  0.01,   day],
    [:e,   0.017,    0.012,  1],
    [:m1,  25.0,     2.3,    m_sun],
    [:K1,  81.4,     1.3,    km_per_s],
]) 

priors = SideKicks.createPriors(
    logm1_dist = Uniform(0.1,3), # in log(Msun)
    logm2_dist = Uniform(0.1,3), # in log(Msun)
    logP_dist  = Uniform(-1,3),   # in log(days)
    vkick_dist = Exponential(1), # in 100 km/s
    frac_dist  = Uniform(0,1.0)
)

##

mcmcStruct = SideKicks.RunKickMCMC(
        pre_supernova_orbit = :circular,
        observations = obs,
        priors = priors,
        nuts_warmup_count = 200,
        nuts_acceptance_rate = 0.8,
        nsamples = 2000,
        nchains = 8)

##

results = mcmcStruct.results

plotting_props_obs_check = SideKicks.createPlottingProps([
    [:m1,    m_sun,    [15,40],        L"M_1\;[M_{\odot}]"],
    [:m2_f,  m_sun,    [6,14],         L"M_{2f}\;[M_{\odot}]"],
    [:P_f,   day,      [10.3,10.5],    L"P_f\;[\mathrm{days}]"],
    [:e_f,   1,        [0,0.1],        L"e_f"],
    [:K1,    km_per_s, [80,85],        L"K_1  \;[\mathrm{km s}^{-1}]"],
])

f = create_corner_plot(results, plotting_props_obs_check,
    tickfontsize=10 ,
    xticklabelrotation=pi/4, 
    show_CIs=true,
    supertitle="Check known quantities"
    )
save("vfts243_obs_check.png", f)

##

plotting_props = SideKicks.createPlottingProps([
    #[:m1,    m_sun,    [15,40],        L"M_1\;[M_{\odot}]"],
    #[:m2_f,  m_sun,    [6,14],         L"M_{2f}\;[M_{\odot}]"],
    #[:P_f,   day,      [10.37,10.43],  L"P_f\;[\mathrm{days}]"],
    #[:e_f,   1,        [0,0.1],        L"e_f"],
    [:m2,    m_sun,     [0,20],        L"M_2  \;[M_{\odot}]"],
    [:dm2,    m_sun,    [0,5],        L"ΔM_2  \;[M_{\odot}]"],
    [:P,     day,      missing,        L"P  \;[\mathrm{days}]"],
    #[:a,     r_sun,    missing,        L"a  \;[R_{odot}]"],
    #[:i_f,   pi,        missing,        L"i_f  \;[π rad]"],
    [:vkick, km_per_s,  [0,50],        L"v_{kick}  \;[\mathrm{km s}^{-1}]"],
    #[:a_f,   r_sun,    missing,        L"a_f  \;[R_{odot}]"],
    #[:K1,    km_per_s, missing,        L"K_1  \;[\mathrm{km s}^{-1}]"],
    #[:K2,    km_per_s, missing,        L"K_2  \;[\mathrm{km s}^{-1}]"],
    #[:frac,  1,        missing,        L"f_{fb}"],
])

f = create_corner_plot(results, plotting_props,
    tickfontsize=10 ,
    xticklabelrotation=pi/4, 
    show_CIs=true,
    supertitle="Derived quantities"
    )
save("vfts243.png", f)