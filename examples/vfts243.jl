using SideKicks
using BenchmarkTools
##
(output_values, loglikelihoods_cauchy, loglikelihoods_normal, weights, chains) = SideKicks.KickMCMC(model_type = :circular,
         observed_properties = [:P, :e, :m1, :K1],
         observed_values = [10.4031*day, 0.017, 25.0*m_sun, 81.4*km_per_s],
         observed_errors = [0.01*day, 0.012, 2.3*m_sun, 1.3*km_per_s],
         nuts_warmup_count = 200,
         nuts_acceptance_rate = 0.8,
         nsamples = 200,
         nchains = 8)
##

#observed_values2 = [10.4031, 0.017, 25.0, 81.4]
#observed_errors2 = [1, 0.012, 2.3, 1.3]
#obs2 = [:P, :e, :m1, :K1]
#
###
#example_model2 = SideKicks.createSimpleCircularMCMCModel(obs2, observed_values2, observed_errors2);
#iterations = 200
###
#@time chain2 = sample(example_model2, NUTS(5_000,0.8), MCMCThreads(), iterations, 8);
#
###
#result = SideKicks.extract_chain(chain2, obs2, observed_values2, observed_errors2, :simple)
#
###
#using CairoMakie
#f = create_corner_plot(result, [:vkick,:ejecta_mass, :m2_f], ["vkick", "ejecta_mass", "BH mass"], [0.68,0.95,0.997], 0.68, Figure(), ranges=[[0,30],[0,3],[5,25]])
#save("vfts243.png", f)
#