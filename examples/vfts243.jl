using SideKicks
using Turing
SideKicks.create_symbolic_functions_list()
##

observed_values2 = [10.4, 0.016, 24.7, 82.84]
observed_errors2 = [0.01, 0.008, 3.9, 0.53]
obs2 = [:P,:e,:m1, :K1]
##

example_model2 = SideKicks.createEccentricMCMCModel(obs2, observed_values2, observed_errors2,
    SideKicks.symbolic_functions_list, e_dist=Uniform(0,1e-5));

iterations = 1000
#chain = sample(example_model2, NUTS(2_000,0.6), MCMCThreads(), iterations, 8);
chain0 = sample(example_model2, NUTS(1_000,0.8), iterations);

##
using BenchmarkTools
function do_a_test(functions_list)
    a_f, e_f, v_N, v_E, v_r, Ω_f, ω_f, ι_f = SideKicks.generalized_post_kick_parameters_a_e(
        1e10,0.1,0.1,0.1,1e33,1e33,0.9*1e33,10*1e7,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,functions_list)
end
@benchmark do_a_test(SideKicks.symbolic_functions_list)

##
example_model2 = SideKicks.createCircularMCMCModel(obs2, observed_values2, observed_errors2,
    SideKicks.symbolic_functions_list);

iterations = 1_000
#chain = sample(example_model2, NUTS(2_000,0.6), MCMCThreads(), iterations, 8);
chain = sample(example_model2, NUTS(1_000,0.8), iterations);

##
example_model2 = SideKicks.createSimpleCircularMCMCModel(obs2, observed_values2, observed_errors2);
iterations = 10000
chain2 = sample(example_model2, NUTS(1_000,0.8), iterations);

##
chain2[:vkick]

using CairoMakie

f = Figure();
ax = Axis(f[1,1])
stephist!(ax, [chain2[:vkick].data[i]*100 for i in 1:10_000], bins=100)
#stephist!(ax, [chain0[:vkick].data[i] for i in 1:1_000 for j in 1:10], bins=70)
f

##
f = Figure();
ax = Axis(f[1,1])
stephist!(ax, [chain2[:frac].data[i] for i in 1:10_000], bins=100)
#stephist!(ax, [chain0[:frac].data[i] for i in 1:1_000 for j in 1:10], bins=70)
f