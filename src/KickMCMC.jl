using Turing
using Distributions

#struct WrappedCauchy <: ContinuousUnivariateDistribution
#    mu::T
#    sigma::Float64
#    WrappedCauchy(mu, sigma) = new(Float64(mu), Float64(sigma))
#end
#
struct WrappedCauchy{T1<:Real, T2<:Real} <: ContinuousUnivariateDistribution
    μ::T1
    σ::T2
end
Distributions.logpdf(d::WrappedCauchy, x::Real) = log(1/(2*π)*(sinh(d.σ)))-log(cosh(d.σ)-cos(x-d.μ))
Distributions.pdf(d::WrappedCauchy, x::Real) = 1/(2*π)*(sinh(d.σ))/(cosh(d.σ)-cos(x-d.μ))

"""
    createEccentricMCMCModel()

Create a Turing model to perform an MCMC sampling of the
pre-explosion and kick properties of a system

#Arguments:
- props: a KickProps variable with necessary prior information

#Output:
- model: A Turing model for sampling
"""
function createEccentricMCMCModel(observations::Vector{Symbol}, observed_values::Vector{Float64}, observed_errors::Vector{Float64},
    functions_list;
    logm1_i_dist::ContinuousUnivariateDistribution = Uniform(0.1,3), 
    logm2_i_dist::ContinuousUnivariateDistribution = Uniform(0.1,3),
    logP_i_dist::ContinuousUnivariateDistribution = Uniform(-1,3),
    e_dist::ContinuousUnivariateDistribution = Uniform(0,0.01),
    vkick_dist::ContinuousUnivariateDistribution = FlatPos(0),
    frac_dist::ContinuousUnivariateDistribution = Uniform(0,1.0))

    # provided observed values
    valid_values = [:P, :e, :K1, :K2, :m1, :m2, :Ω, :ω, :ι, :v_N, :v_E, :v_r]
    for obs ∈ observations
        if obs ∉ valid_values
            throw(DomainError(obs, "Allowed observations are only [:P, :e, :K1, :K2, :m1, :m2, :Ω, :ω, :ι, :v_N, :v_E, :v_r]"))
        end
    end

    @model function kick_model(obs::Vector{Symbol}, obs_vals::Vector{Float64}, obs_errs::Vector{Float64})
        # set priors
        #Pre-explosion masses and orbital period
        logm1_i ~ logm1_i_dist
        m1_i = 10^(logm1_i)
        logm2_i ~ logm2_i_dist
        m2_i = 10^(logm2_i)
        logP_i ~ logP_i_dist
        P_i = 10^(logP_i)
        a_i = kepler_a_from_P(P_i,m1_i,m2_i)
        e ~ e_dist
        cosι ~ Uniform(0,1)
        sinι = sqrt(1-cosι^2)
        #azimuthal angles are computed by sampling random points with a circularly symmetric distribution
        xν ~ Normal()
        yν ~ Normal()
        normν = 1/sqrt(xν^2+yν^2)
        cosν = xν*normν
        sinν = yν*normν
        xΩ ~ Normal()
        yΩ ~ Normal()
        normΩ = 1/sqrt(xΩ^2+yΩ^2)
        cosΩ = xΩ*normΩ
        sinΩ = yΩ*normΩ
        xω ~ Normal()
        yω ~ Normal()
        normω = 1/sqrt(xω^2+yω^2)
        cosω = xω*normω
        sinω = yω*normω

        #Post-explosion masses
        frac ~ frac_dist
        m2_f = m2_i*frac # star 2 explodes, star 1 is kept fixed

        #Kick parameters
        vkick ~ vkick_dist
        cosθ ~ Uniform(-1,1)
        sinθ = sqrt(1-cosθ^2)
        xϕ ~ Normal(0,1)
        yϕ ~ Normal(0,1)
        normϕ = 1/sqrt(xϕ^2+yϕ^2)
        cosϕ = xϕ*normϕ
        sinϕ = yϕ*normϕ

        #m1 is assumed to remain constant
        a_f, e_f, v_N, v_E, v_r, Ω_f, ω_f, ι_f = 
            generalized_post_kick_parameters_a_e(
                        a_i,e,sinν,cosν,m1_i*m_sun,m2_i*m_sun,m2_f*m_sun,vkick*1e7,sinθ,cosθ,sinϕ,cosϕ,sinΩ,cosΩ,sinω,cosω,sinι,cosι,functions_list)
        P_f = kepler_P_from_a(a_f,m1_i,m2_f)
        
        #@show "in model", (P_f,e_f)
        K1 = RV_semiamplitude_K(P_f, e_f, ι_f, m1_i, m2_f)
        K2 = RV_semiamplitude_K(P_f, e_f, ι_f, m2_f, m1_i)

        for (i, obs_symbol) in enumerate(obs)
            if obs_symbol == :P
                obs_vals[i] ~ Cauchy(P_f, obs_errs[i])
            elseif obs_symbol == :e
                obs_vals[i] ~ Cauchy(e_f, obs_errs[i])
            elseif obs_symbol == :K1
                obs_vals[i] ~ Cauchy(K1, obs_errs[i])
            elseif obs_symbol == :K2
                obs_vals[i] ~ Cauchy(K2, obs_errs[i])
            elseif obs_symbol == :m1
                obs_vals[i] ~ Cauchy(m1_i, obs_errs[i])
            elseif obs_symbol == :m2
                obs_vals[i] ~ Cauchy(m2_f, obs_errs[i])
            elseif obs_symbol == :Ω
                #obs_vals[i] ~ VonMises(Ω_f, 1/obs_errs[i]^2)#Normal(Ω_f, obs_errs[i])
                obs_vals[i] ~ WrappedCauchy(Ω_f, obs_errs[i])
            elseif obs_symbol == :ω
                #obs_vals[i] ~ VonMises(ω_f, 1/obs_errs[i]^2)#Normal(ω_f, obs_errs[i])
                obs_vals[i] ~ WrappedCauchy(ω_f, obs_errs[i])
            elseif obs_symbol == :ι
                obs_vals[i] ~ Cauchy(ι_f, obs_errs[i])
            elseif obs_symbol == :v_N
                obs_vals[i] ~ Cauchy(v_N, obs_errs[i])
            elseif obs_symbol == :v_E
                obs_vals[i] ~ Cauchy(v_E, obs_errs[i])
            elseif obs_symbol == :v_r
                obs_vals[i] ~ Cauchy(v_r, obs_errs[i])
            end
        end

        return m1_i, m2_i, m2_f, P_i, vkick, P_f, e_f, K1, K2, v_N, v_E, v_r, Ω_f, ω_f, ι_f
    end

    return kick_model(observations, observed_values, observed_errors)

end
