using Turing
using Distributions
using StatsPlots

# The distribution of the true velocity, this is what we want to infer
σ_v = 10.0
μ_v = 35.0
v_dist = Normal(μ_v, σ_v)

# we assume the measurement errors themselves follow an independent distribution
σ_e = 1
μ_e = 5
v_error_dist = Normal(μ_e, σ_e)

# we make n_sample "measurements", we want to infer the distribution of the true velocity with this
n_sample = 10_000
true_velocities = rand(v_dist,n_sample)
errors = abs.(rand(v_error_dist, n_sample))
measured_velocities = zeros(n_sample)
for i in 1:n_sample
    measured_v_error_dist = Normal(0.0, errors[i])
    measured_velocities[i] = true_velocities[i] + rand(measured_v_error_dist)
end



##

# Create the model and sample it. This assumes the measured velocity is the sum of two things,
# the real velocity and the error. Likelihood is determined considering that the distribution
# of the sum of two random variables is another normally distributed random variable.
@model function v_model(measured_velocities, errors)
    σ_v ~ Uniform(0.0,200.0)
    μ_v ~ Uniform(0.0,200.0)

    for i in eachindex(measured_velocities)
        new_σ = sqrt(σ_v^2 + errors[i]^2)
        measured_velocities[i] ~ Normal(μ_v, new_σ)
    end
end
iterations = 1_000
model = v_model(measured_velocities, errors)
chain = sample(model, NUTS(10_000,0.8), MCMCThreads(), iterations, 8);

##
#f=Figure()

#Axis(f, chain)
#scatter!(chain)
#plot(chain)
chain
