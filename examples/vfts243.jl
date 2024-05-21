using SideKicks
using BenchmarkTools
using Logging

obs = SideKicks.Observations(
    props        =  [:P,       :e,     :m1,       :K1],                           
    vals         =  [10.4031,  0.017,  25.0,      81.4],                          
    errs         =  [0.01,     0.012,  2.3,       1.3],                           
    units        =  [day,      1,      m_sun,     km_per_s],
    #names_latex  =  ["Per",    "Ecc",  "Mass_1",  "K1"],                           
    #units_latex  =  ["d",      "-",    "M_odot",  "km_per_s"]  #  Come  back  to  this
) 

println(obs.props)

##

#Logging.disable_logging(Logging.Info)

@benchmark mcmcStruct = SideKicks.RunKickMCMC(
        pre_supernova_orbit = :circular,
        observations = obs,
        nuts_warmup_count = 200,
        nuts_acceptance_rate = 0.8,
        nsamples = 200,
        nchains = 8)

##

# test
results = mcmcStruct.results
results[:e_f]
size(results[:e_f])

println(results[:P][1:10])
println(results[:vkick][1:10])
println(results[:m2_f][1:10])
##

using CairoMakie
plotting_props = SideKicks.PlottingProps(
    props  = [:m1,   :m2_f, :P_f, :e_f],
    units = [m_sun, m_sun, day, 1, ],
    ranges = [[0, 25], [0, .0005], [0, 1], [0, 1]], 
    names_latex = ["m1",   "m2_f", "P_f", "e_f" ]
)

#=
plotting_props = SideKicks.PlottingProps(
    props = [:vkick, :P, :m2_f], 
    units = [km_per_s, day, m_sun],
    ranges = missing, #[[0,100],[0,15],[5,25]], 
    names_latex = ["vkick km/s", "P_0", "M2"], # RTW clean up
#    props = [:vkick, :ejecta_mass, :m2_f], 
#    units = [km_per_s, m_sun, m_sun],
#    ranges = [[0,30],[0,3],[5,25]], 
#    names_latex = ["vkick km/s", "M_ej", "M2"], # RTW clean up
)
=#

#f = create_corner_plot(results, [:vkick, :ejecta_mass, :m2_f], ["vkick", "ejecta_mass", "BH mass"], [0.68,0.95,0.997], 0.68, Figure(), ranges=[[0,30],[0,3],[5,25]])
f = create_corner_plot(results, plotting_props)
save("vfts243.png", f)


