using Turing

"""
    createMCMCModel(props::KickProps)

Create a Turing model to perform an MCMC sampling of the
pre-explosion and kick properties of a system

#Arguments:
- props: a KickProps variable with necessary prior information

#Output:
- model: A Turing model for sampling
"""
function createMCMCModel(observations::Vector{Symbol}, observed_values::Vector{Float64}, observed_errors::Vector{Float64};
    have_i::Bool = false,
    vkick_div_vrel_dist::ContinuousUnivariateDistribution = Uniform(0,1+sqrt(2)), 
    logm1_i_dist::ContinuousUnivariateDistribution = Uniform(0.1,3), 
    logm2_i_dist::ContinuousUnivariateDistribution = Uniform(0.1,3),
    frac_dist::ContinuousUnivariateDistribution = Uniform(0,1.0),
    cos_θ_dist::ContinuousUnivariateDistribution = Uniform(-1,1), 
    ϕ_dist::ContinuousUnivariateDistribution = Uniform(0,2*π),
    cos_i_dist::ContinuousUnivariateDistribution = Uniform(0,1),
    i_dist::ContinuousUnivariateDistribution = Uniform(0,π/2),
    logP_i_dist::ContinuousUnivariateDistribution = Uniform(-1,3))

    # provided observed values
    valid_values = [:P, :e, :K1, :K2, :vsys]
    for obs ∈ observations
        if obs ∉ valid_values
            throw(DomainError(obs, "Allowed observations are only [:P, :e, :K1, :K2, :vsys]"))
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

        #Post-explosion masses
        frac ~ frac_dist
        m2_f = m2_i*frac # star 2 explodes, star 1 is kept fixed
        mtilde = (m1_i+m2_f)/(m1_i+m2_i)
        vrel = sqrt(cgrav*(m1_i+m2_i)*m_sun/a_i)

        #Kick parameters
        vkick_div_vrel ~ vkick_div_vrel_dist
        vkick = vkick_div_vrel*vrel
        cos_θ ~ cos_θ_dist
        θ = acos(cos_θ)
        ϕ ~ ϕ_dist

        #Post-explosion orbital inclination
        if (have_i)
            i ~ i_dist
        else
            cos_i ~ cos_i_dist
            i = acos(cos_i)
        end

        #m1 is assumed to remain constant
        P_f, e_f = post_kick_parameters_P(P_i,mtilde,vkick_div_vrel,θ,ϕ)
        vsys = systemic_velocity(vrel,m1_i,m2_i,m2_f,vkick,θ,ϕ)
        #@show "in model", (P_f,e_f)
        K1 = RV_semiamplitude_K(P_f, e_f, i, m1_i, m2_f)
        K2 = RV_semiamplitude_K(P_f, e_f, i, m2_f, m1_i)

        for (i, obs_symbol) in enumerate(obs)
            if obs_symbol == :P
                obs_vals[i] ~ Normal(P_f, obs_errs[i])
            elseif obs_symbol == :e
                obs_vals[i] ~ Normal(e_f, obs_errs[i])
            elseif obs_symbol == :K1
                obs_vals[i] ~ Normal(K1, obs_errs[i])
            elseif obs_symbol == :K2
                obs_vals[i] ~ Normal(K2, obs_errs[i])
            else
                obs_vals[i] ~ Normal(vsys, obs_errs[i])

            end
        end

        return m1_i, m2_i, m2_f, P_i, vkick_div_vrel*vrel,θ,ϕ, P_f, e_f, K1, K2, vsys
    end

    return kick_model(observations, observed_values, observed_errors)

end
