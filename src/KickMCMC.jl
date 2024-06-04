using Turing
using Distributions

########################
### 
### Struct definitions
### 
########################

"""
    struct WrappedCauchy{T1<:Real, T2<:Real} <: ContinuousUnivariateDistribution

The WrappedCauchy distribution resembles the Cauchy distribution defined on the unit 
circle from 0 to 2π, with the endpoints wrapped back to each other.
"""
struct WrappedCauchy{T1<:Real, T2<:Real} <: ContinuousUnivariateDistribution
    μ::T1
    σ::T2
end
Distributions.logpdf(d::WrappedCauchy, x::Real) = log(1/(2*π)*(sinh(d.σ)))-log(cosh(d.σ)-cos(x-d.μ))
Distributions.pdf(d::WrappedCauchy, x::Real) = 1/(2*π)*(sinh(d.σ))/(cosh(d.σ)-cos(x-d.μ))

"""
    mutable struct Observations

Observations contains the symbols, values, errors, and units of each observed parameter.
"""
@kwdef mutable struct Observations
    props::Vector{Symbol}
    vals::Vector{Float64}
    errs::Vector{Float64}
    units::Vector{Float64}
end

"""
    mutable struct Priors

Priors contains the prior distribution of each of the desired parameters

RTW: check that N and E are identical to RA and Dec...
"""
@kwdef mutable struct Priors
    logm1_dist::ContinuousUnivariateDistribution
    logm2_dist::ContinuousUnivariateDistribution
    logP_dist::ContinuousUnivariateDistribution 
    vkick_dist::ContinuousUnivariateDistribution
    frac_dist::ContinuousUnivariateDistribution 
    e_dist::ContinuousUnivariateDistribution       
    rv_env_dist::ContinuousUnivariateDistribution 
    pmra_env_dist::ContinuousUnivariateDistribution  
    pmdec_env_dist::ContinuousUnivariateDistribution 
    parallax_dist::ContinuousUnivariateDistribution 
    # RTW make this parallax
end

"""
    mutable struct KickMCMCResults

KickMCMCResults contains the MCMC model and Observations structs, the Results dict,
the chains that resulted from the MCMC, and the parameters that went into the sampler.
"""
@kwdef mutable struct KickMCMCResults
    # TODO RTW: do I need types for everything? what do I use if non-trivial?
    mcmc_model # RTW todo
    observations::Observations
    results::Dict{Symbol, Matrix{Float64}}
    chains # RTW todo
    nuts_warmup_count::Int
    nuts_acceptance_rate::Float64
    nsamples::Int
end





########################
### 
### Function definitions
### 
########################

function createObservations(obs_matrix::Vector{Vector{Any}})
    obs_matrix = stack(obs_matrix) # make into actual matrix
    return Observations(
        props = obs_matrix[1,:],
        vals  = obs_matrix[2,:],
        errs  = obs_matrix[3,:],
        units = obs_matrix[4,:])
end

function createPriors(;
    logm1_dist = Uniform(0.1,3), # in log(Msun)
    logm2_dist = Uniform(0.1,3), # in log(Msun)
    logP_dist  = Uniform(-1,3),   # in log(days)
    vkick_dist = Exponential(1), # in 100 km/s
    frac_dist  = Uniform(0,1.0))
    # For ecc systems - ignored for circular case
    e_dist = Uniform(0,0.01)
    # RTW: fix these priors, they are too narrow for what we know about these quantities...
    rv_env_dist    = Normal(0,10) # in 100 km/s
    pmra_env_dist  = Normal(0,10) # in mas/yr 
    pmdec_env_dist = Normal(0,10) # in mas/yr 
    parallax_dist  = Normal(0.01,0.005) # in 1/kpc
    return Priors(
        logm1_dist = logm1_dist,  
        logm2_dist = logm2_dist,  
        logP_dist  = logP_dist ,  
        vkick_dist = vkick_dist,  
        frac_dist  = frac_dist,
        e_dist = e_dist,
        rv_env_dist    = rv_env_dist,
        pmra_env_dist  = pmra_env_dist,
        pmdec_env_dist = pmdec_env_dist,
        parallax_dist  = parallax_dist)
end

"""
    createSimpleCircularMCMCModel(observations, observed_values, observed_errors)

Description
Create a Turing model to perform an MCMC sampling of the pre-explosion 
and kick properties of a system, assuming pre-explosion circularity.

RTW this is more simplistic than just using a circular model, it's also
assuming you know the eccentricity and don't care about radial velocity etc.

# Arguments:
- observations:    the parameters taken from observations [Vector{Symbol}]
- observed_values: the values of the parameters           [Vector{Float64}] 
- observed_errors: the errors of the observations         [Vector{Float64}]

# Output:
- kickmodel: A Turing model for sampling
"""
function createCircularMCMCModel(;
    observations::Observations,
    priors::Priors,
    likelihood = :Cauchy,
    bhModel = arbitraryEjectaBH
    )

    logm1_dist = priors.logm1_dist
    logm2_dist = priors.logm2_dist
    logP_dist = priors.logP_dist
    vkick_dist = priors.vkick_dist
    frac_dist = priors.frac_dist
    rv_env_dist = priors.rv_env_dist
    pmra_env_dist = priors.pmra_env_dist
    pmdec_env_dist = priors.pmdec_env_dist
    parallax_dist = priors.parallax_dist

    # RTW: why do we need this?
    ## provided observed values
    #valid_values = [:P, :e, :K1, :K2, :m1, :m2]
    #for obs ∈ observations.props
    #    if obs ∉ valid_values
    #        throw(DomainError(obs, "Allowed observations are only [:P, :e, :K1, :K2, :m1, :m2]"))
    #    end
    #end
    if !(likelihood == :Cauchy || likelihood == :Normal)
        throw(DomainError(likelihood, "likelihood must be either :Cauchy or :Normal"))
    end

    @model function create_MCMC_model(obs_vals, obs_errs) 

        # set priors
        #Pre-explosion masses and orbital period
        logm1 ~ logm1_dist
        m1 = 10^(logm1)*m_sun
        logm2 ~ logm2_dist
        m2 = 10^(logm2)*m_sun
        logP ~ logP_dist
        P = 10^(logP)*day
        a = kepler_a_from_P(m1=m1, m2=m2, P=P)
        cosi ~ Uniform(0,1)
        i_f = acos(cosi)

        #Post-explosion masses
        frac ~ frac_dist
        m2_f = bhModel(m2, frac) # star 2 explodes, star 1 is kept fixed

        #Kick parameters
        vkick_100kms ~ vkick_dist  
        vkick = vkick_100kms *100*km_per_s
        cosθ ~ Uniform(-1,1)
        θ = acos(cosθ)
        xϕ ~ Normal(0,1)
        yϕ ~ Normal(0,1)
        normϕ = 1/sqrt(xϕ^2 + yϕ^2)
        cosϕ = xϕ*normϕ
        ϕ = acos(cosϕ)

        #m1 is assumed to remain constant
        a_f, e_f = post_supernova_circular_orbit_a(m1=m1, m2=m2, a=a, m2_f=m2_f, vkick=vkick, θ=θ, ϕ=ϕ)
        P_f = kepler_P_from_a(m1=m1, m2=m2_f, a=a_f)
        K1 = RV_semiamplitude_K1(m1=m1, m2=m2_f, P=P_f, e=e_f, i=i_f)
        K2 = RV_semiamplitude_K1(m1=m2_f, m2=m1, P=P_f, e=e_f, i=i_f)

        likelihood == :Cauchy ?
            likelihood_dist = Cauchy :
            likelihood_dist = Normal
        obs_err = obs_errs[i]

        for i in eachindex(observations.props)
            obs_symbol = observations.props[i]
            if obs_symbol == :P
                param = P_f
            elseif obs_symbol == :e
                param = e_f
            elseif obs_symbol == :K1
                param = K1
            elseif obs_symbol == :K2
                param = K2
            elseif obs_symbol == :m1
                param = m1
            elseif obs_symbol == :m2
                param = m2_f
            else
                continue
            end

            obs_vals[i] ~ likelihood_dist(param, obs_err) # Needs to have this [i] in the obs_vals, no idea why...
        end

        # other params
        dm2 = m2 - m2_f
        return     (m1,    m2,    P,   a,     i_f, vkick, m2_f,  a_f,   P_f,  e_f,  K1,       K2, frac, dm2)
    end
    return_props = [:m1,   :m2,   :P,  :a,    :i_f, :vkick, :m2_f, :a_f,  :P_f, :e_f, :K1,      :K2, :frac, :dm2]

    obs_vals_cgs = observations.vals .* observations.units
    obs_errs_cgs = observations.errs .* observations.units
    return [create_MCMC_model(obs_vals_cgs, obs_errs_cgs), return_props]
end


"""
    createEccentricMCMCModel(observations, observed_values, observed_errors)

Create a Turing model to perform an MCMC sampling of the pre-explosion 
and kick properties of a system, assuming pre-explosion eccentricity.

# RTW does arccos just work? Do I need to worry about domain/range issues?
# Check velocities, the conversions are a bit funky

# Arguments:
- observations:    the parameters taken from observations [Vector{Symbol}]
- observed_values: the values of the parameters           [Vector{Float64}] 
- observed_errors: the errors of the observations         [Vector{Float64}]

# Output:
- kickmodel: A Turing model for sampling
"""
function createEccentricMCMCModel(;
    observations::Observations,
    priors::Priors,
    likelihood = :Cauchy,
    bhModel = arbitraryEjectaBH)
    
    logm1_dist = priors.logm1_dist
    logm2_dist = priors.logm2_dist
    logP_dist = priors.logP_dist
    e_dist = priors.e_dist
    vkick_dist = priors.vkick_dist
    frac_dist = priors.frac_dist
    rv_env_dist = priors.rv_env_dist
    pmra_env_dist = priors.pmra_env_dist
    pmdec_env_dist = priors.pmdec_env_dist
    parallax_dist = priors.parallax_dist

    # provided observed values
    # RTW: why do we need this?
    #valid_values = [:P, :e, :K1, :K2, :m1, :m2, :Ω, :ω, :i, :v_N, :v_E, :v_r]
    #for obs ∈ observations.props
    #    if obs ∉ valid_values
    #        throw(DomainError(obs, "Allowed observations are only [:P, :e, :K1, :K2, :m1, :m2, :Ω, :ω, :i, :v_N, :v_E, :v_r]"))
    #    end
    #end
    if !(likelihood == :Cauchy || likelihood == :Normal)
        throw(DomainError(likelihood, "likelihood must be either :Cauchy or :Normal"))
    end

    @model function create_MCMC_model(obs_props, obs_vals, obs_errs) 
        # set priors
        #Pre-explosion masses and orbital period
        logm1 ~ logm1_dist
        m1 = 10^(logm1)*m_sun
        logm2 ~ logm2_dist
        m2 = 10^(logm2)*m_sun
        logP ~ logP_dist
        P = 10^(logP)*days
        a = kepler_a_from_P(m1=m1, m2=m2, P=P)
        e ~ e_dist
        cosi ~ Uniform(0,1)
        sini = sqrt(1-cosi^2)
        #azimuthal angles are computed by sampling random points with a circularly symmetric distribution
        #we take all true anomalies to be equally likely. This is corrected by weighting later
        xν ~ Normal()
        yν ~ Normal()
        normν = 1/sqrt(xν^2+yν^2)
        cosν = xν*normν
        ν = arccos(cosν)
        xΩ ~ Normal()
        yΩ ~ Normal()
        normΩ = 1/sqrt(xΩ^2+yΩ^2)
        cosΩ = xΩ*normΩ
        Ω = arccos(cosΩ)
        xω ~ Normal()
        yω ~ Normal()
        normω = 1/sqrt(xω^2+yω^2)
        cosω = xω*normω
        ω = arccos(cosω)

        #rv_env_dist = 
        #pmra_env_dist =
        #pmdec_env_dist 
        #parallax_dist =
        # TODO: how to specify either pm/parallax or v_Sys priors?


        #Post-explosion masses
        frac ~ frac_dist
        m2_f = bhModel(m2, frac) # star 2 explodes, star 1 is kept fixed

        #Kick parameters
        vkick_100kms ~ vkick_dist
        vkick = vkick_100kms*100*km_per_s
        cosθ ~ Uniform(-1,1)
        θ = arccos(cosθ)
        xϕ ~ Normal(0,1)
        yϕ ~ Normal(0,1)
        normϕ = 1/sqrt(xϕ^2+yϕ^2)
        cosϕ = xϕ*normϕ
        ϕ = arccos(cosϕ)
        if yϕ < 0
            ϕ = 2π - ϕ
        end

        #Initial systemic velocity parameters
        v_N_100kms ~ vsys_N_dist
        v_E_100kms ~ vsys_E_dist
        v_r_100kms ~ vsys_r_dist
        v_N = v_N_100kms*100*km_per_s
        v_E = v_E_100kms*100*km_per_s
        v_r = v_r_100kms*100*km_per_s

        #m1 is assumed to remain constant, no impact velocity
        a_f, e_f, Ω_f, ω_f, i_f, v_N, v_E, v_r = 
            post_supernova_general_orbit_parameters(m1=m1, m2=m2, a=a, e=e, m2_f=m2_f, 
                vkick=vkick, θ=θ, ϕ=ϕ, ν=ν, Ω=Ω, ω=ω, i=i)
        P_f = kepler_P_from_a(m1=m1, m2=m2_f, a=a_f)
        K1 = RV_semiamplitude_K1(m1=m1, m2=m2_f, P=P_f, e=e_f, i=i_f)
        K2 = RV_semiamplitude_K1(m1=m2_f, m2=m1, P=P_f, e=e_f, i=i_f)

        likelihood == :Cauchy ?
            likelihood_dist = Cauchy :
            likelihood_dist = Normal
        obs_err = obs_errs[i]

        # RTW TODO: when to use VonMises vs WrappedCauchy?
        for i in eachindex(obs_props)
            obs_symbol = obs_props[i]
            if obs_symbol == :P
                param = P_f
            elseif obs_symbol == :e
                param = e_f
            elseif obs_symbol == :K1
                param = K1
            elseif obs_symbol == :K2
                param = K2
            elseif obs_symbol == :m1
                param = m1
            elseif obs_symbol == :m2
                param = m2_f
            elseif obs_symbol == :i
                param = i_f
            elseif obs_symbol == :Ω
                param = Ω_f
                if likelihood == :Cauchy 
                    likelihood_dist = WrappedCauchy
                else
                    likelihood_dist = VonMises
                    obs_err = 1/obs_errs[i]^2
                end
            elseif obs_symbol == :ω
                param = ω_f
                if likelihood == :Cauchy 
                    likelihood_dist = WrappedCauchy
                else
                    likelihood_dist = VonMises
                    obs_err = 1/obs_errs[i]^2
                end
            # RTW: what is going on with these, v_N - v_N?
            #elseif obs_symbol == :v_N
            #    param = v_N - v_N
            #elseif obs_symbol == :v_E
            #    param = v_E - v_E
            #elseif obs_symbol == :v_r
            #    param = v_r - v_r
            else
                continue
            end

            obs_vals[i] ~ likelihood_dist(param, obs_err)
        end

        return     (m1,    m2,    P,   a,     i_f, vkick, m2_f,  a_f,   P_f,  e_f,  K1,       K2, frac)
    end
    # RTW : fix the props later
    return_props = [:m1,   :m2,   :P,  :a,    :i_f, :vkick, :m2_f, :a_f,  :P_f, :e_f, :K1,      :K2, :frac]

    # Need to combine some of the observations to compare against the predicted output
    obs_vals_cgs = observations.vals .* observations.units
    obs_errs_cgs = observations.errs .* observations.units




    return [create_MCMC_model(obs_vals_cgs, obs_errs_cgs), return_props]
    
end


function RunKickMCMC(; pre_supernova_orbit, observations::Observations, priors::Priors,
        nuts_warmup_count, nuts_acceptance_rate, nsamples, nchains)

    if (pre_supernova_orbit==:circular)
        mcmc_cauchy, props_cauchy = SideKicks.createCircularMCMCModel( observations=observations, priors=priors, likelihood=:Cauchy)
        mcmc_normal, props_normal = SideKicks.createCircularMCMCModel( observations=observations, priors=priors, likelihood=:Normal)
    elseif (pre_supernova_orbit==:eccentric)
        mcmc_cauchy, props_cauchy = SideKicks.createEccentricMCMCModel( observations=observations, priors=priors, likelihood=:Cauchy)
        mcmc_normal, props_normal = SideKicks.createEccentricMCMCModel( observations=observations, priors=priors, likelihood=:Normal)
    else
        throw(ArgumentError("pre_supernova_orbit=:$pre_supernova_orbit is an invalid option. Can be either :circular or :eccentric"))
    end

    #run the MCMC - this is the slow step!
    chains = sample(mcmc_cauchy,
                    NUTS(nuts_warmup_count,nuts_acceptance_rate),
                    MCMCThreads(),
                    nsamples,
                    nchains);

    # compute weights to get sampling from a normal distribution
    loglikelihoods_cauchy = pointwise_loglikelihoods(mcmc_cauchy, chains)
    loglikelihoods_normal = pointwise_loglikelihoods(mcmc_normal, chains)

    # Obtain the generated values from the chains
    output_values = generated_quantities(mcmc_cauchy, chains)
    
    # Combine into a dictionary of matrices, where the keys are props and the matrices are (chains x samples)
    results = Dict() 
    for i_prop in eachindex(props_cauchy)                                                    
        prop = props_cauchy[i_prop]
        chain_array = zeros(Float64, nchains, nsamples) # Matrix (nchain x nsample)
        for i_chain in 1:nchains
            for i_sample in 1:nsamples
                chain_array[i_chain, i_sample] = output_values[i_sample, i_chain][i_prop]
            end
        end
        results[prop] = chain_array
    end
    # Add weights to dict
    logweights = zeros(Float64, nchains, nsamples) # Matrix (nchain x nsample)
    for i_chain in 1:nchains
        for i_sample in 1:nsamples
            for dict_key in keys(loglikelihoods_cauchy)
                logweights[i_chain, i_sample] = loglikelihoods_normal[dict_key][i_sample, i_chain]
                                               -loglikelihoods_cauchy[dict_key][i_sample, i_chain]
            end
        end
    end
    logweights = logweights .- maximum(logweights) # set max weights = 1
    # RTW something broke here
    weights = exp.(logweights)
    if (all(isfinite(weights)))
        results[:weights] = weights 
    else
        println("Weights are off!!")
        results[:weights] = ones(size(weights)) 
        # TODO: throw an exception here
    end

    #RTW this was from the previous iteration, is this still relevant?
    # if we did a general model, we need to weight the true anomaly
    #if model_type==:general
    #    res[:weight] .= res[:weight].*sqrt.(1 .- res[:e_f].^2).^3 ./ (1 .+ res[:e_f].*cos.(res[:ν])).^2
    #end
   
    return KickMCMCResults(
        mcmc_model = mcmc_cauchy, 
        observations = observations,
        results = results,
        chains = chains,
        nuts_warmup_count = nuts_warmup_count,
        nuts_acceptance_rate = nuts_acceptance_rate,
        nsamples = nsamples
    )
end


# for the 1D distributions in the corner plots, plot each chain (maybe optionally)
# 2D distributions might get a bit messy with this

