using Turing
using Distributions

"""
RTW TODO: figure out correct way to document structs
"""
struct WrappedCauchy{T1<:Real, T2<:Real} <: ContinuousUnivariateDistribution
    μ::T1
    σ::T2
end
Distributions.logpdf(d::WrappedCauchy, x::Real) = log(1/(2*π)*(sinh(d.σ)))-log(cosh(d.σ)-cos(x-d.μ))
Distributions.pdf(d::WrappedCauchy, x::Real) = 1/(2*π)*(sinh(d.σ))/(cosh(d.σ)-cos(x-d.μ))


"""
    createSimpleCircularMCMCModel(observations, observed_values, observed_errors)

Description
Create a Turing model to perform an MCMC sampling of the pre-explosion 
and kick properties of a system, assuming pre-explosion circularity.

# Arguments:
- observations:    the parameters taken from observations [Vector{Symbol}]
- observed_values: the values of the parameters           [Vector{Float64}] 
- observed_errors: the errors of the observations         [Vector{Float64}]

# Output:
- kickmodel: A Turing model for sampling
"""
function createCircularMCMCModel(observations::Vector{Symbol}, observed_values::Vector{Float64}, observed_errors::Vector{Float64};
    bhModel = arbitraryEjectaBH,
    logm1_dist::ContinuousUnivariateDistribution = Uniform(0.1,3), # in log(Msun)
    logm2_dist::ContinuousUnivariateDistribution = Uniform(0.1,3), # in log(Msun)
    logP_dist::ContinuousUnivariateDistribution = Uniform(-1,3),   # in log(days)
    vkick_dist::ContinuousUnivariateDistribution = Exponential(1), # in 100 km/s
    frac_dist::ContinuousUnivariateDistribution = Uniform(0,1.0),
    likelihood = :Cauchy)

    # provided observed values
    valid_values = [:P, :e, :K1, :K2, :m1, :m2]
    for obs ∈ observations
        if obs ∉ valid_values
            throw(DomainError(obs, "Allowed observations are only [:P, :e, :K1, :K2, :m1, :m2]"))
        end
    end
    if !(likelihood == :Cauchy || likelihood == :Normal)
        throw(DomainError(obs, "likelihood must be either :Cauchy or :Normal"))
    end

    @model function kick_model(obs::Vector{Symbol}, obs_vals::Vector{Float64}, obs_errs::Vector{Float64})
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

        for i in eachindex(obs)
            obs_symbol = obs[i]
            if obs_symbol == :P
                likelihood == :Cauchy ?
                    obs_vals[i] ~ Cauchy(P_f, obs_errs[i]) :
                    obs_vals[i] ~ Normal(P_f, obs_errs[i])
            elseif obs_symbol == :e
                likelihood == :Cauchy ?
                    obs_vals[i] ~ Cauchy(e_f, obs_errs[i]) :
                    obs_vals[i] ~ Normal(e_f, obs_errs[i])
            elseif obs_symbol == :K1
                likelihood == :Cauchy ?
                    obs_vals[i] ~ Cauchy(K1, obs_errs[i]) :
                    obs_vals[i] ~ Normal(K1, obs_errs[i])
            elseif obs_symbol == :K2
                likelihood == :Cauchy ?
                    obs_vals[i] ~ Cauchy(K2, obs_errs[i]) :
                    obs_vals[i] ~ Normal(K2, obs_errs[i])
            elseif obs_symbol == :m1
                likelihood == :Cauchy ?
                    obs_vals[i] ~ Cauchy(m1, obs_errs[i]) :
                    obs_vals[i] ~ Normal(m1, obs_errs[i])
            elseif obs_symbol == :m2
                likelihood == :Cauchy ?
                    obs_vals[i] ~ Cauchy(m2_f, obs_errs[i]) :
                    obs_vals[i] ~ Normal(m2_f, obs_errs[i])
            end
        end
        return (m1, m2, P, a, i_f, m2_f, a_f, P_f, e_f, K1, K2)
    end

    return (kick_model(observations, observed_values, observed_errors),
           [:m1, :m2, :P, :a, :i_f, :m2_f, :a_f, :P_f, :e_f, :K1, :K2]) 

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
function createEccentricMCMCModel(observations::Vector{Symbol}, observed_values::Vector{Float64}, observed_errors::Vector{Float64};
    bhModel = arbitraryEjectaBH,
    logm1_dist::ContinuousUnivariateDistribution = Uniform(0.1,3), # in log(Msun)
    logm2_dist::ContinuousUnivariateDistribution = Uniform(0.1,3), # in log(Msun)
    logP_dist::ContinuousUnivariateDistribution = Uniform(-1,3),   # in log(days)
    e_dist::ContinuousUnivariateDistribution = Uniform(0,0.01),
    vkick_dist::ContinuousUnivariateDistribution = Exponential(1), # in 100 km/s
    vsys_N_dist::ContinuousUnivariateDistribution = Normal(0,0.1), # in 100 km/s
    vsys_E_dist::ContinuousUnivariateDistribution = Normal(0,0.1), # in 100 km/s
    vsys_r_dist::ContinuousUnivariateDistribution = Normal(0,0.1), # in 100 km/s
    frac_dist::ContinuousUnivariateDistribution = Uniform(0,1.0),
    likelihood = :Cauchy)

    # provided observed values
    valid_values = [:P, :e, :K1, :K2, :m1, :m2, :Ω, :ω, :i, :v_N, :v_E, :v_r]
    for obs ∈ observations
        if obs ∉ valid_values
            throw(DomainError(obs, "Allowed observations are only [:P, :e, :K1, :K2, :m1, :m2, :Ω, :ω, :i, :v_N, :v_E, :v_r]"))
        end
    end

    @model function kick_model(obs::Vector{Symbol}, obs_vals::Vector{Float64}, obs_errs::Vector{Float64})
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

        #Post-explosion masses
        frac ~ frac_dist
        m2_f = bhModel(m2, frac) # star 2 explodes, star 1 is kept fixed

        #Kick parameters
        vkick ~ vkick_dist*100*km_per_s
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
        v_N ~ vsys_N_dist*100*km_per_s 
        v_E ~ vsys_E_dist*100*km_per_s 
        v_r ~ vsys_r_dist*100*km_per_s 

        #m1 is assumed to remain constant, no impact velocity
        a_f, e_f, Ω_f, ω_f, i_f, v_N, v_E, v_r = 
            post_supernova_general_orbit_parameters(m1=m1, m2=m2, a=a, e=e, m2_f=m2_f, 
                vkick=vkick, θ=θ, ϕ=ϕ, ν=ν, Ω=Ω, ω=ω, i=i)
        P_f = kepler_P_from_a(m1=m1, m2=m2_f, a=a_f)
        K1 = RV_semiamplitude_K1(m1=m1, m2=m2_f, P=P_f, e=e_f, i=i_f)
        K2 = RV_semiamplitude_K1(m1=m2_f, m2=m1, P=P_f, e=e_f, i=i_f)

        for i in eachindex(obs)
            obs_symbol = obs[i]
            if obs_symbol == :P
                obs_vals[i] ~ likelihood_dist(P_f, obs_errs[i])
            elseif obs_symbol == :e
                obs_vals[i] ~ likelihood_dist(e_f, obs_errs[i])
            elseif obs_symbol == :K1
                obs_vals[i] ~ likelihood_dist(K1, obs_errs[i])
            elseif obs_symbol == :K2
                obs_vals[i] ~ likelihood_dist(K2, obs_errs[i])
            elseif obs_symbol == :m1
                obs_vals[i] ~ likelihood_dist(m1, obs_errs[i])
            elseif obs_symbol == :m2
                obs_vals[i] ~ likelihood_dist(m2_f, obs_errs[i])
            elseif obs_symbol == :Ω
                #obs_vals[i] ~ VonMises(Ω_f, 1/obs_errs[i]^2)#Normal(Ω_f, obs_errs[i])
                obs_vals[i] ~ angular_likelihood_dist(Ω_f, obs_errs[i])
            elseif obs_symbol == :ω
                #obs_vals[i] ~ VonMises(ω_f, 1/obs_errs[i]^2)#Normal(ω_f, obs_errs[i])
                obs_vals[i] ~ angular_likelihood_dist(ω_f, obs_errs[i])
            elseif obs_symbol == :i
                obs_vals[i] ~ likelihood_dist(i_f, obs_errs[i])
            elseif obs_symbol == :v_N
                obs_vals[i] ~ likelihood_dist(v_N - v_N, obs_errs[i])
            elseif obs_symbol == :v_E
                obs_vals[i] ~ likelihood_dist(v_E - v_E, obs_errs[i])
            elseif obs_symbol == :v_r
                obs_vals[i] ~ likelihood_dist(v_r - v_r, obs_errs[i])
            end
        end
    end
    
    return kick_model(observations, observed_values, observed_errors)
    
end
    
mutable struct KickMCMC
    model_type::Symbol
    chains::Chains
    result::Vector{Dict{Symbol, Vector{Float64}}}
    nuts_warmup_count::Int
    nuts_acceptance_rate::Float64
    nsamples::Int
    observed_properties::Vector{Symbol}
    observed_values::Vector{Float64}
    observed_errors::Vector{Float64}
end

function KickMCMC(;model_type, observed_properties, observed_values, observed_errors,
    nuts_warmup_count, nuts_acceptance_rate, nsamples, nchains)
    if !(model_type==:circular || model_type==:eccentric)
        throw(ArgumentError("model_type=:$model_type is an invalid option. Can be either :circular or :eccentric"))
    end
    if (model_type==:circular)
        (model_cauchy, output_names) = SideKicks.createCircularMCMCModel(
            observed_properties, observed_values, observed_errors, likelihood=:Cauchy);
        (model_normal, output_names) = SideKicks.createCircularMCMCModel(
            observed_properties, observed_values, observed_errors, likelihood=:Normal);
    else
        (model_cauchy, output_names) = SideKicks.createEccentricMCMCModel(
            observed_properties, observed_values, observed_errors, likelihood=:Cauchy);
        (model_normal, output_names) = SideKicks.createEccentricMCMCModel(
            observed_properties, observed_values, observed_errors, likelihood=:Normal);
    end

    #run the MCMC
    chains = sample(model_cauchy,
                    NUTS(nuts_warmup_count,nuts_acceptance_rate),
                    MCMCThreads(),
                    nsamples,
                    nchains);

    # compute weights to get sampling from a normal distribution
    loglikelihoods_cauchy = pointwise_loglikelihoods(model_cauchy, chains)
    loglikelihoods_normal = pointwise_loglikelihoods(model_normal, chains)
    weights = [zeros(nsamples) for i in 1:nchains]
    for i_chain in 1:nchains
        for i_sample in 1:nsamples
            for dict_key in keys(loglikelihoods_cauchy)
                weights[i_chain][i_sample] -= loglikelihoods_cauchy[dict_key][i_sample, i_chain]
                weights[i_chain][i_sample] += loglikelihoods_normal[dict_key][i_sample, i_chain]
            end
            weights[i_chain][i_sample] = exp(weights[i_chain][i_sample])
        end
    end

    # obtain the generated values from the chains
    output_values = generated_quantities(model_cauchy, chains)
    
    return (output_values, loglikelihoods_cauchy, loglikelihoods_normal, weights, chains)
end
