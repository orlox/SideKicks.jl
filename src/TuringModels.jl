using Turing
using Distributions

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
    struct ModVonMises{T1<:Real, T2<:Real} <: ContinuousUnivariateDistribution

This is just a wrapper on top of the VonMises distribution (as defined in Distributions.jl)
to extend its domain. This is because the domain of VonMises is defined to be
[μ-π, μ+π], and the angles we are concerned with range from [0,2π]
"""
struct ModVonMises{T1<:Real, T2<:Real, T3<:ContinuousUnivariateDistribution} <: ContinuousUnivariateDistribution
    μ::T1
    κ::T2
    vonMisesDist::T3
end
function ModVonMises(μ, κ)
    return ModVonMises(μ, κ, VonMises(μ, κ))
end
function Distributions.logpdf(d::ModVonMises, x::Real)
    # If x is outside the range [μ-π, μ+π], we need to shift it by the correct
    # amount of 2π to fit it there
    if x > d.μ + π
        return logpdf(d.vonMisesDist, x-2π)
    elseif x < d.μ -π
        return logpdf(d.vonMisesDist, x+2π)
    else
        return logpdf(d.vonMisesDist,x)
    end
end
function Distributions.pdf(d::ModVonMises, x::Real)
    # If x is outside the range [μ-π, μ+π], we need to shift it by the correct
    # amount of 2π to fit it there
    if x > d.μ + π
        return pdf(d.vonMisesDist, x-2π)
    elseif x < d.μ -π
        return pdf(d.vonMisesDist, x+2π)
    else
        return pdf(d.vonMisesDist,x)
    end
end

"""
    create_simplified_mcmc_model(observations, observed_values, observed_errors)

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
function create_simplified_mcmc_model(;
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

    valid_values = [:P_f, :e_f, :K1, :K2, :m1_f, :m2_f, :i_f]
    for prop ∈ observations.props
        if prop ∉ valid_values
            throw(DomainError(observation.props, "Allowed observations are only [:P_f, :e_f, :K1, :K2, :m1_f, :m2_f, :i_f]"))
        end
    end

    if !(likelihood == :Cauchy || likelihood == :Normal)
        throw(DomainError(likelihood, "likelihood must be either :Cauchy or :Normal"))
    end

    @model function create_mcmc_model(props, obs_vals, obs_errs) 

        # set priors
        #Pre-explosion masses and orbital period
        logm1 ~ logm1_dist
        m1_i = 10^(logm1)*m_sun
        logm2 ~ logm2_dist
        m2_i = 10^(logm2)*m_sun
        logP ~ logP_dist
        P_i = 10^(logP)*day
        a_i = kepler_a_from_P(m1=m1_i, m2=m2_i, P=P_i)
        cosi ~ Uniform(0,1)
        i_f = acos(cosi)

        #Post-explosion masses
        frac ~ frac_dist
        m2_f = bhModel(m2_i, frac) # star 2 explodes, star 1 is kept fixed

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
        m1_f = m1_i
        a_f, e_f = post_supernova_circular_orbit_a(m1_i=m1_i, m2_i=m2_i, a_i=a_i, m2_f=m2_f, vkick=vkick, θ=θ, ϕ=ϕ)
        P_f = kepler_P_from_a(m1=m1_f, m2=m2_f, a=a_f)
        K1 = RV_semiamplitude_K1(m1=m1_f, m2=m2_f, P=P_f, e=e_f, i=i_f)
        K2 = RV_semiamplitude_K1(m1=m2_f, m2=m1_f, P=P_f, e=e_f, i=i_f)

        use_cauchy = likelihood == :Cauchy
        for ii in eachindex(props)
            obs_symbol = props[ii]
            if obs_symbol == :P_f
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(P_f, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(P_f, obs_errs[ii])
            elseif obs_symbol == :e_f
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(e_f, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(e_f, obs_errs[ii])
            elseif obs_symbol == :K1
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(K1, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(K1, obs_errs[ii])
            elseif obs_symbol == :K2
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(K2, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(K2, obs_errs[ii])
            elseif obs_symbol == :m1_f
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(m1_f, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(m1_f, obs_errs[ii])
            elseif obs_symbol == :m2_f
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(m2_f, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(m2_f, obs_errs[ii])
            elseif obs_symbol == :i_f
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(i_f, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(i_f, obs_errs[ii])
            end
        end

        # other params
        dm2 = m2_i - m2_f
        return     (m1_i,    m2_i,    P_i,   a_i,     i_f,  vkick,  m2_f,  a_f,   P_f,  e_f,  K1,  K2, frac, dm2) 
    end
    return_props = [:m1_i,   :m2_i,   :P_i,  :a_i,    :i_f, :vkick, :m2_f, :a_f,  :P_f, :e_f, :K1, :K2, :frac, :dm2] 

    obs_vals_cgs = observations.vals .* observations.units
    obs_errs_cgs = observations.errs .* observations.units
    return [create_mcmc_model(observations.props, obs_vals_cgs, obs_errs_cgs), return_props]
end


"""
    create_general_mcmc_model(observations, observed_values, observed_errors)

Create a Turing model to perform an MCMC sampling of the pre-explosion 
and kick properties of a system, assuming pre-explosion eccentricity.

# RTW does acos just work? Do I need to worry about domain/range issues?
# Check velocities, the conversions are a bit funky

# Arguments:
- observations:    the parameters taken from observations [Vector{Symbol}]
- observed_values: the values of the parameters           [Vector{Float64}] 
- observed_errors: the errors of the observations         [Vector{Float64}]

# Output:
- kickmodel: A Turing model for sampling
"""
function create_general_mcmc_model(;
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
    venv_N_100kms_dist = priors.venv_N_100kms_dist
    venv_E_100kms_dist = priors.venv_E_100kms_dist
    venv_r_100kms_dist = priors.venv_r_100kms_dist

    valid_values = [:P_f, :e_f, :K1, :K2, :m1_f, :m2_f, :Ω_f, :ω_f, :i_f, :v_N, :v_E, :v_r]
    for prop ∈ observations.props
        if prop ∉ valid_values
            throw(DomainError(observation.props, "Allowed observations are only [:P_f, :e_f, :K1, :K2, :m1_f, :m2_f, :Ω_f, :ω_f, :i_f, :v_N, :v_E, :v_r]"))
        end
    end
    if !(likelihood == :Cauchy || likelihood == :Normal)
        throw(DomainError(likelihood, "likelihood must be either :Cauchy or :Normal"))
    end

    @model function create_mcmc_model(props, obs_vals, obs_errs) 
        # set priors
        #Pre-explosion masses and orbital period
        logm1 ~ logm1_dist
        m1_i = 10^(logm1)*m_sun
        logm2 ~ logm2_dist
        m2_i = 10^(logm2)*m_sun
        logP ~ logP_dist
        P_i = 10^(logP)*day
        a_i = kepler_a_from_P(m1=m1_i, m2=m2_i, P=P_i)
        e_i ~ e_dist

        #azimuthal angles are computed by sampling random points with a circularly symmetric distribution
        #we take all true anomalies to be equally likely. This is corrected by weighting later
        xν ~ Normal()
        yν ~ Normal()
        normν = 1/sqrt(xν^2+yν^2)
        cosν = xν*normν
        ν_i = acos(cosν)
        xΩ ~ Normal()
        yΩ ~ Normal()
        normΩ = 1/sqrt(xΩ^2+yΩ^2)
        cosΩ = xΩ*normΩ
        Ω_i = acos(cosΩ)
        if yΩ < 0
            Ω_i = 2π - Ω_i
        end
        xω ~ Normal()
        yω ~ Normal()
        normω = 1/sqrt(xω^2+yω^2)
        cosω = xω*normω
        ω_i = acos(cosω)
        if yω < 0
            ω_i = 2π - ω_i
        end
        cosi ~ Uniform(0,1)
        i_i = acos(cosi)

        #Post-explosion masses
        frac ~ frac_dist
        m2_f = bhModel(m2_i, frac) # star 2 explodes, star 1 is kept fixed

        #Kick parameters
        vkick_100kms ~ vkick_dist
        vkick = vkick_100kms*100*km_per_s
        cosθ ~ Uniform(-1,1)
        θ = acos(cosθ)
        xϕ ~ Normal(0,1)
        yϕ ~ Normal(0,1)
        normϕ = 1/sqrt(xϕ^2+yϕ^2)
        cosϕ = xϕ*normϕ
        ϕ = acos(cosϕ)
        if yϕ < 0
            ϕ = 2π - ϕ
        end

        #Initial systemic velocity parameters
        venv_N_100kms ~ venv_N_100kms_dist
        venv_E_100kms ~ venv_E_100kms_dist
        venv_r_100kms ~ venv_r_100kms_dist
        venv_N = venv_N_100kms*100*km_per_s 
        venv_E = venv_E_100kms*100*km_per_s 
        venv_r = venv_r_100kms*100*km_per_s 

        #m1 is assumed to remain constant, no impact velocity - TODO: relax this later
        m1_f = m1_i
        a_f, e_f, Ω_f, ω_f, i_f, vsys_N, vsys_E, vsys_r = 
            post_supernova_general_orbit_parameters(m1_i=m1_i, m2_i=m2_i, a_i=a_i, e_i=e_i, m2_f=m2_f, 
                vkick=vkick, θ=θ, ϕ=ϕ, ν_i=ν_i, Ω_i=Ω_i, ω_i=ω_i, i_i=i_i)
        P_f = kepler_P_from_a(m1=m1_f, m2=m2_f, a=a_f)
        K1 = RV_semiamplitude_K1(m1=m1_f, m2=m2_f, P=P_f, e=e_f, i=i_f)
        K2 = RV_semiamplitude_K1(m1=m2_f, m2=m1_f, P=P_f, e=e_f, i=i_f)
        vsys = sqrt( vsys_N^2 + vsys_E^2 + vsys_r^2)

        v_N = venv_N + vsys_N
        v_E = venv_E + vsys_E
        v_r = venv_r + vsys_r
        
        use_cauchy = likelihood == :Cauchy
        for ii in eachindex(props)
            obs_symbol = props[ii]
            if obs_symbol == :P_f
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(P_f, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(P_f, obs_errs[ii])
            elseif obs_symbol == :e_f
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(e_f, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(e_f, obs_errs[ii])
            elseif obs_symbol == :K1
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(K1, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(K1, obs_errs[ii])
            elseif obs_symbol == :K2
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(K2, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(K2, obs_errs[ii])
            elseif obs_symbol == :m1_f
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(m1_f, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(m1_f, obs_errs[ii])
            elseif obs_symbol == :m2_f
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(m2_f, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(m2_f, obs_errs[ii])
            elseif obs_symbol == :i_f
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(i_f, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(i_f, obs_errs[ii])
            elseif obs_symbol == :Ω_f
                use_cauchy ?
                    obs_vals[ii] ~ WrappedCauchy(Ω_f, obs_errs[ii]) :
                    obs_vals[ii] ~ ModVonMises(Ω_f, 1/obs_errs[ii]^2)
            elseif obs_symbol == :ω_f
                use_cauchy ?
                    obs_vals[ii] ~ WrappedCauchy(ω_f, obs_errs[ii]) :
                    obs_vals[ii] ~ ModVonMises(Ω_f, 1/obs_errs[ii]^2)
            elseif obs_symbol == :v_N
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(v_N, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(v_N, obs_errs[ii])
            elseif obs_symbol == :v_E
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(v_E, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(v_E, obs_errs[ii])
            elseif obs_symbol == :v_r
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(v_r, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(v_r, obs_errs[ii])
            end
        end

        dm2 = m2_i - m2_f
        return     ( m1_i,  m2_i,  P_i,  e_i,  a_i,  ν_i,  vkick,  frac,  dm2,  i_f,  ω_f,  Ω_f,  m1_f, m2_f,  a_f,  P_f,  e_f,  K1,  K2,  v_N,  v_E,  v_r,  vsys)
    end
    return_props = [:m1_i, :m2_i, :P_i, :e_i, :a_i, :ν_i, :vkick, :frac, :dm2, :i_f, :ω_f, :Ω_f, :m1_f,:m2_f, :a_f, :P_f, :e_f, :K1, :K2, :v_N, :v_E, :v_r, :vsys]

    # Need to combine some of the observations to compare against the predicted output
    obs_vals_cgs = observations.vals .* observations.units
    obs_errs_cgs = observations.errs .* observations.units

    return [create_mcmc_model(observations.props, obs_vals_cgs, obs_errs_cgs), return_props]
    
end
