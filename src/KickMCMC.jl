using Turing
using Distributions
using Logging

"""
    mutable struct KickMCMC

KickMCMC contains the MCMC model and Observations structs, the Results dict,
the chains that resulted from the MCMC, and the parameters that went into the sampler.
"""
mutable struct KickMCMC
    mcmc_model_cauchy 
    mcmc_model_normal 
    observations::Observations
    observations_string::String
    priors::Priors
    priors_string::String
    results::Dict{Symbol, Matrix{Float64}}
    chains 
    nuts_warmup_count::Int
    nuts_acceptance_rate::Float64
    nsamples::Int
end

function KickMCMC(; which_model, observations::Tuple{Observations, String}, priors::Tuple{Priors, String},
        nuts_warmup_count, nuts_acceptance_rate, nsamples, nchains)

    if (which_model==:simplified)
        mcmc_cauchy, props_cauchy = SideKicks.create_simplified_mcmc_model( observations=observations[1], priors=priors[1], likelihood=:Cauchy)
        mcmc_normal, props_normal = SideKicks.create_simplified_mcmc_model( observations=observations[1], priors=priors[1], likelihood=:Normal)
    elseif (which_model==:general)
        mcmc_cauchy, props_cauchy = SideKicks.create_general_mcmc_model( observations=observations[1], priors=priors[1], likelihood=:Cauchy)
        mcmc_normal, props_normal = SideKicks.create_general_mcmc_model( observations=observations[1], priors=priors[1], likelihood=:Normal)
    else
        throw(ArgumentError("which_model=:$which_model is an invalid option. Can be either :simplified or :general"))
    end

    # Run the MCMC - this is the slow step!
    chains = sample(mcmc_cauchy,
                    NUTS(nuts_warmup_count,nuts_acceptance_rate),
                    MCMCThreads(),
                    nsamples,
                    nchains);

    # pointwise_loglikelihood can give a lot of spurious warnings. We suppress them
    # by doing this part with a specific logger.
    (loglikelihoods_cauchy, loglikelihoods_normal) = with_logger(ConsoleLogger(Error)) do
        # Compute weights to get sampling from a normal distribution
        loglikelihoods_cauchy = pointwise_loglikelihoods(mcmc_cauchy, chains)
        loglikelihoods_normal = pointwise_loglikelihoods(mcmc_normal, chains)
        return (loglikelihoods_cauchy, loglikelihoods_normal)
    end

    # Obtain the generated values from the chains
    output_values = generated_quantities(mcmc_cauchy, chains)
    
    # Combine into a dictionary of matrices, where the keys are props and the matrices are (chains x samples)
    results = Dict() 
    for i_prop in eachindex(props_cauchy)                                                    
        prop = props_cauchy[i_prop]
        chain_array = zeros(Float64, nsamples, nchains) # Matrix (nsample x nchain)
        for i_chain in 1:nchains
            for i_sample in 1:nsamples
                chain_array[i_sample, i_chain] .= output_values[i_sample, i_chain][i_prop]
            end
        end
        results[prop] = chain_array
    end
    # Add weights to dict
    logweights = zeros(Float64, nchains, nsamples) # Matrix (nsample x nchain)
    for i_chain in 1:nchains
        for i_sample in 1:nsamples
            for dict_key in keys(loglikelihoods_cauchy)
                logweights[i_sample, i_chain] += 
                   loglikelihoods_normal[dict_key][i_chain, i_chain] - loglikelihoods_cauchy[dict_key][i_chain, i_chain]
            end
        end
    end
    logweights = logweights .- maximum(logweights) # set max weights = 1
    weights = exp.(logweights)
    if (all(isfinite(weights)))
        results[:weights] = weights 
    else
        throw(ErrorException("Failed to reweight samples"))
    end

    # for the general model, need to reweight by the true anomaly
    if which_model==:general
        results[:weights] .= results[:weights].*sqrt.(1 .- results[:e_i].^2).^3 ./ (1 .+ results[:e_i].*cos.(results[:ν_i])).^2
    end
    
    return KickMCMC(
        mcmc_cauchy, 
        mcmc_normal, 
        observations[1],
        observations[2],
        priors[1],
        priors[2],
        results,
        chains,
        nuts_warmup_count,
        nuts_acceptance_rate,
        nsamples
    )
end
