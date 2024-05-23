using SideKicks
using BenchmarkTools
using CairoMakie
using Distributions

obs = SideKicks.Observations(
    props        =  [:P,       :e,     :m1,       :K1],                           
    vals         =  [10.4031,  0.017,  25.0,      81.4],                          
    errs         =  [0.01   ,  0.012,  2.3 ,       1.3],                           
    units        =  [day,      1,      m_sun,     km_per_s],
) 

priors = SideKicks.createPriors(
    logm1_dist = Uniform(0.1,3), # in log(Msun)
    logm2_dist = Uniform(0.1,3), # in log(Msun)
    logP_dist  = Uniform(-1,3),   # in log(days)
    vkick_dist = Exponential(1), # in 100 km/s
    frac_dist  = Uniform(0,1.0)
)

println(obs.props)
##

mcmcStruct = SideKicks.RunKickMCMC(
        pre_supernova_orbit = :circular,
        observations = obs,
        priors = priors,
        nuts_warmup_count = 200,
        nuts_acceptance_rate = 0.8,
        nsamples = 200,
        nchains = 8)

##

results = mcmcStruct.results

plotting_props = SideKicks.PlottingProps(
    props  = [:m1,   :m2_f, :P_f, :e_f],
    units = [m_sun, m_sun, day, 1 ],
    ranges = [[15, 45], [6, 14], [9, 11], [0, 1]], 
    names_latex = ["m1",   "m2_f", "P_f", "e_f" ]
)

f = create_corner_plot(results, plotting_props)
save("vfts243.png", f)


