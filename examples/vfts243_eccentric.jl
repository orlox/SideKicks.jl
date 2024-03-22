using SideKicks
using Turing
using BenchmarkTools
##

observed_values2 = [10.4031, 0.017, 25.0, 81.4]
observed_errors2 = [0.01, 0.012, 2.3, 1.3]
obs2 = [:P,:e,:m1, :K1]

##
example_model2 = SideKicks.createEccentricMCMCModel(obs2, observed_values2, observed_errors2;e_dist = Uniform(0,0.1));
iterations = 20_000
@time chain2 = sample(example_model2, NUTS(5_000,0.8), MCMCThreads(), iterations, 8);
#@time chain2 = sample(example_model2, NUTS(5_000,0.8), iterations);

##
result = SideKicks.extract_chain(chain2, obs2, observed_values2, observed_errors2, :general)

##
using CairoMakie
f = create_corner_plot(result, [:vkick,:ejecta_mass, :m2_f], ["vkick", "ejecta_mass", "BH mass"], [0.68,0.95,0.997], 0.68, Figure(), ranges=[[0,30],[0,3],[5,25]])
save("vfts243.png", f)