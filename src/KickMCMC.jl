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
    logm1_i_dist::ContinuousUnivariateDistribution = Uniform(0.1,3), 
    logm2_i_dist::ContinuousUnivariateDistribution = Uniform(0.1,3),
    logP_i_dist::ContinuousUnivariateDistribution = Uniform(-1,3),
    e_dist::ContinuousUnivariateDistribution = Uniform(0,0.01),
    vkick_dist::ContinuousUnivariateDistribution = Exponential(1), # in 100 km/s
    vsys_N_i_dist::ContinuousUnivariateDistribution = Normal(0,0.1), # in 100 km/s
    vsys_E_i_dist::ContinuousUnivariateDistribution = Normal(0,0.1), # in 100 km/s
    vsys_r_i_dist::ContinuousUnivariateDistribution = Normal(0,0.1), # in 100 km/s
    frac_dist::ContinuousUnivariateDistribution = Uniform(0,1.0))

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
        logm1_i ~ logm1_i_dist
        m1_i = 10^(logm1_i)
        logm2_i ~ logm2_i_dist
        m2_i = 10^(logm2_i)
        logP_i ~ logP_i_dist
        P_i = 10^(logP_i)
        a_i = kepler_a_from_P(;m_1=m1_i, m_2=m2_i, P=P_i)
        e_i ~ e_dist
        cosi ~ Uniform(0,1)
        sini = sqrt(1-cosi^2)
        #azimuthal angles are computed by sampling random points with a circularly symmetric distribution
        #we take all true anomalies to be equally likely. This is corrected by weighting later
        xν ~ Normal()
        yν ~ Normal()
        normν = 1/sqrt(xν^2+yν^2)
        cosν_i = xν*normν
        ν_i = arccos(cosν_i)
        #sinν = yν*normν
        xΩ ~ Normal()
        yΩ ~ Normal()
        normΩ = 1/sqrt(xΩ^2+yΩ^2)
        cosΩ_i = xΩ*normΩ
        Ω_i = arccos(cosΩ_i)
        #sinΩ = yΩ*normΩ
        xω ~ Normal()
        yω ~ Normal()
        normω = 1/sqrt(xω^2+yω^2)
        cosω_i = xω*normω
        ω_i = arccos(cosω_i)
        #sinω = yω*normω

        #Post-explosion masses
        frac ~ frac_dist
        m2_f = bhModel(m2_i, frac) # star 2 explodes, star 1 is kept fixed

        #Kick parameters
        vkick ~ vkick_dist
        # RTW what units is this in?
        cosθ ~ Uniform(-1,1)
        θ = arccos(θ)
        #sinθ = sqrt(1-cosθ^2)
        xϕ ~ Normal(0,1)
        yϕ ~ Normal(0,1)
        normϕ = 1/sqrt(xϕ^2+yϕ^2)
        cosϕ = xϕ*normϕ
        ϕ = arccos(ϕ)
        #sinϕ = yϕ*normϕ

        #Initial systemic velocity parameters
        v_N_i ~ vsys_N_i_dist*1e7         # in cm/s
        v_E_i ~ vsys_E_i_dist*1e7         # in cm/s
        v_r_i ~ vsys_r_i_dist*1e7         # in cm/s

        #m1 is assumed to remain constant, no impact velocity
        a_f, e_f, Ω_f, ω_f, i_f, v_n, v_e, v_rad = 
            post_supernova_general_orbit_parameters(;m_1i=m1_i, m_2i=m2_i, a_i=a_i, e_i=e_i, m_2f=m2_f, 
                vkick=vkick*100, θ=θ, ϕ=ϕ, ν_i=ν_i, Ω_i=Ω_i, ω_i=ω_i, i_i=i_i)
        P_f = kepler_P_from_a(;m_1=m1_i, m_2=m2_f, a=a_f)
        K1 = RV_semiamplitude_K1(;m_1=m1_i, m_2=m2_f, P=P_f, e=e_f, i=i_f)
        K2 = RV_semiamplitude_K1(;m_1=m2_f, m_2=m1_i, P=P_f, e=e_f, i=i_f)

        for i in eachindex(obs)
            obs_symbol = obs[i]
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
            elseif obs_symbol == :i
                obs_vals[i] ~ Cauchy(i_f, obs_errs[i])
            elseif obs_symbol == :v_N
                obs_vals[i] ~ Cauchy(v_N - v_N_i, obs_errs[i])
            elseif obs_symbol == :v_E
                obs_vals[i] ~ Cauchy(v_E - v_E_i, obs_errs[i])
            elseif obs_symbol == :v_r
                obs_vals[i] ~ Cauchy(v_r - v_r_i, obs_errs[i])
            end
        end
    end
    
    return kick_model(observations, observed_values, observed_errors)
    
end
    
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
function createSimpleCircularMCMCModel(observations::Vector{Symbol}, observed_values::Vector{Float64}, observed_errors::Vector{Float64};
    bhModel = arbitraryEjectaBH,
    logm1_i_dist::ContinuousUnivariateDistribution = Uniform(0.1,3), 
    logm2_i_dist::ContinuousUnivariateDistribution = Uniform(0.1,3),
    logP_i_dist::ContinuousUnivariateDistribution = Uniform(-1,3),
    vkick_dist::ContinuousUnivariateDistribution = Exponential(1),
    frac_dist::ContinuousUnivariateDistribution = Uniform(0,1.0))

    # provided observed values
    valid_values = [:P, :e, :K1, :K2, :m1, :m2]
    for obs ∈ observations
        if obs ∉ valid_values
            throw(DomainError(obs, "Allowed observations are only [:P, :e, :K1, :K2, :m1, :m2]"))
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
        a_i = kepler_a_from_P(;m_1=m1_i, m_2=m2_i, P=P_i)
        cosi ~ Uniform(0,1)
        i_f = acos(cosi)

        #Post-explosion masses
        frac ~ frac_dist
        m2_f = bhModel(m2_i, frac) # star 2 explodes, star 1 is kept fixed

        #Kick parameters
        vkick ~ vkick_dist
        cosθ ~ Uniform(-1,1)
        θ = acos(cosθ)
        xϕ ~ Normal(0,1)
        yϕ ~ Normal(0,1)
        normϕ = 1/sqrt(xϕ^2+yϕ^2)
        cosϕ = xϕ*normϕ
        ϕ = acos(cosϕ)

        #m1 is assumed to remain constant
        a_f, e_f = post_supernova_circular_orbit_a(;m_1i=m1_i, m_2i=m2_i, a_i=a_i, m_2f=m2_f, vkick=100*vkick, θ=θ, ϕ=ϕ)
        P_f = kepler_P_from_a(;m_1=m1_i, m_2=m2_i, a=a_i)
        K1 = RV_semiamplitude_K1(;m_1=m1_i, m_2=m2_f, P=P_f, e=e_f, i=i_f)

        for i in eachindex(obs)
            obs_symbol = obs[i]
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
            end
        end
    end

    return kick_model(observations, observed_values, observed_errors)

end

"""
    extract_chain(chain, observations, observed_values, observed_errors,
        model_type; bhModel = arbitraryEjectaBH)

#TODO Description

# Arguments:
#TODO
- chain:   
- observations:   
- observed_values:   
- observed_errors:   

# Output:
#TODO
- res:
"""
function extract_chain(chain, observations, observed_values, observed_errors,
        model_type; bhModel = arbitraryEjectaBH)
    #simple function to change ranges from [-π,π] to [0,2π]
    shift_range = x-> x<0 ? x+2*π : x


    res = Dict()
    res[:logm1_i] = reduce(vcat,chain[:logm1_i])
    res[:m1_i] = 10 .^ res[:logm1_i]
    res[:logm2_i] = reduce(vcat,chain[:logm2_i])
    res[:m2_i] = 10 .^ res[:logm2_i]
    res[:logP_i] = reduce(vcat,chain[:logP_i])
    res[:P_i] = 10 .^ res[:logP_i]
    res[:a_i] = kepler_a_from_P.(;m_1=res[:m1_i], m_2=res[:m2_i], P=res[:P_i])
    if model_type==:simple
        res[:e_i] = zeros(length(res[:a_i]))
    elseif model_type==:general
        res[:e_i] = reduce(vcat,chain[:e])
    else
        throw(ArgumentError("model_type=:$model_type is an invalid option"))
    end
    if model_type==:general
        res[:ν_i] = reduce(vcat,shift_range.([atan(y,x) for (x,y) in zip(chain[:xν],chain[:yν])]))
        res[:Ω_i] = reduce(vcat,shift_range.([atan(y,x) for (x,y) in zip(chain[:xΩ],chain[:yΩ])]))
        res[:ω_i] = reduce(vcat,shift_range.([atan(y,x) for (x,y) in zip(chain[:xω],chain[:yω])]))
        res[:i_i] = reduce(vcat,[acos(x) for x in chain[:cosi]])
    end
    res[:frac] = reduce(vcat,chain[:frac])
    res[:m2_f] = bhModel.(res[:m2_i],res[:frac])
    res[:ejecta_mass] = res[:m2_i] .- res[:m2_f]
    res[:vkick] = reduce(vcat,chain[:vkick])*100 # MCMC is done in units of 100 km/s, turn into km/s
    res[:θ] = reduce(vcat,[acos(x) for x in chain[:cosθ]])
    res[:ϕ] = reduce(vcat,shift_range.([atan(y,x) for (x,y) in zip(chain[:xϕ],chain[:yϕ])]))
    if model_type==:general
        res[:v_N_i] = reduce(vcat,chain[:v_N_i])
        res[:v_E_i] = reduce(vcat,chain[:v_E_i])
        res[:v_r_i] = reduce(vcat,chain[:v_r_i])
    end 

    sample_num = length(res[:a_i])
    if model_type==:general
        res[:a_f] = Vector{Float64}(undef, sample_num)
        res[:e_f] = Vector{Float64}(undef, sample_num)
        res[:v_N] = Vector{Float64}(undef, sample_num)
        res[:v_E] = Vector{Float64}(undef, sample_num)
        res[:v_r] = Vector{Float64}(undef, sample_num)
        res[:Ω_f] = Vector{Float64}(undef, sample_num)
        res[:ω_f] = Vector{Float64}(undef, sample_num)
        res[:i_f] = Vector{Float64}(undef, sample_num)
        for i in 1:sample_num
            a_f, e_f, Ω_f, ω_f, i_f, v_n, v_e, v_rad = 
            post_supernova_general_orbit_parameters(; m_1i= res[:m1_i][i], m_2i= res[:m2_i][i],
                a_i= res[:a_i][i], e_i= res[:e_i][i], m_2f= res[:m2_f][i], vkick= res[:vkick][i], 
                θ= res[:θ][i], ϕ= res[:ϕ][i], ν_i= res[:ν_i][i], Ω_i= res[:Ω_i][i], ω_i= res[:ω_i][i], i_i= res[:i_i][i])

            res[:a_f][i] = a_f
            res[:e_f][i] = e_f
            res[:v_N][i] = v_N + res[:v_N_i][i]
            res[:v_E][i] = v_E + res[:v_E_i][i]
            res[:v_r][i] = v_r + res[:v_r_i][i]
            res[:Ω_f][i] = Ω_f
            res[:ω_f][i] = ω_f
            res[:i_f][i] = i_f
        end
    elseif model_type==:simple
        res[:a_f] = Vector{Float64}(undef, sample_num)
        res[:e_f] = Vector{Float64}(undef, sample_num)
        for i in 1:sample_num
            a_f, e_f = post_supernova_circular_orbit_a(; m_1i=res[:m1_i][i], m_2i=res[:m2_i][i], 
                                                         a_i=res[:a_i][i], m_2f=res[:m2_f][i], vkick=res[:vkick][i], θ=res[:θ][i], ϕ=res[:ϕ][i])
            res[:a_f][i] = a_f
            res[:e_f][i] = e_f
        end
        res[:i_f] = reduce(vcat,[acos(x) for x in chain[:cosi]])
    end
    res[:P_f] = kepler_P_from_a.(; m_1=res[:m1_i], m_2=res[:m2_f], a=res[:a_f])
    res[:K1] = RV_semiamplitude_K1.(;m_1=res[:m1_i], m_2=res[:m2_f], P=res[:P_f], e=res[:e_f], i=res[:i_f])
    res[:K2] = RV_semiamplitude_K1.(;m_1=res[:m2_f], m_2=res[:m1_i], P=res[:P_f], e=res[:e_f], i=res[:i_f])
    #compute weights in log first
    res[:weight] = zeros(Float64,sample_num)
    for j in 1:sample_num
        for (i, obs_symbol) in enumerate(observations)
            if obs_symbol == :P
                dist1 = Cauchy(res[:P_f][j], observed_errors[i])
                dist2 = Normal(res[:P_f][j], observed_errors[i])
                res[:weight][j] += logpdf(dist2,observed_values[i])-logpdf(dist1,observed_values[i])
            elseif obs_symbol == :e
                dist1 = Cauchy(res[:e_f][j], observed_errors[i])
                dist2 = Normal(res[:e_f][j], observed_errors[i])
                res[:weight][j] += logpdf(dist2,observed_values[i])-logpdf(dist1,observed_values[i])
            elseif obs_symbol == :K1
                dist1 = Cauchy(res[:K1][j], observed_errors[i])
                dist2 = Normal(res[:K1][j], observed_errors[i])
                res[:weight][j] += logpdf(dist2,observed_values[i])-logpdf(dist1,observed_values[i])
            elseif obs_symbol == :K2
                dist1 = Cauchy(res[:K2][j], observed_errors[i])
                dist2 = Normal(res[:K2][j], observed_errors[i])
                res[:weight][j] += logpdf(dist2,observed_values[i])-logpdf(dist1,observed_values[i])
            elseif obs_symbol == :m1
                dist1 = Cauchy(res[:m1_i][j], observed_errors[i])
                dist2 = Normal(res[:m1_i][j], observed_errors[i])
                res[:weight][j] += logpdf(dist2,observed_values[i])-logpdf(dist1,observed_values[i])
            elseif obs_symbol == :m2
                dist1 = Cauchy(res[:m2_f][j], observed_errors[i])
                dist2 = Normal(res[:m2_f][j], observed_errors[i])
                res[:weight][j] += logpdf(dist2,observed_values[i])-logpdf(dist1,observed_values[i])
            elseif obs_symbol == :Ω
                dist1 = WrappedCauchy(res[:Ω_f][j], observed_errors[i])
                dist2 = VonMises(res[:Ω_f][j], 1/observed_errors[i]^2)
                # VonMises works only between μ-π and μ+π, so need to adjust accordingly
                diff = observed_values[i]-res[:Ω_f][j]
                if diff<-π
                    diff += 2*π
                elseif diff>π
                    diff -= 2*π
                end
                res[:weight][j] += logpdf(dist2,res[:Ω_f][j]+diff)-logpdf(dist1,observed_values[i])
            elseif obs_symbol == :ω
                dist1 = WrappedCauchy(res[:ω_f][j], observed_errors[i])
                dist2 = VonMises(res[:ω_f][j], 1/observed_errors[i]^2)
                # VonMises works only between μ-π and μ+π, so need to adjust accordingly
                diff = observed_values[i]-res[:ω_f][j]
                if diff<-π
                    diff += 2*π
                elseif diff>π
                    diff -= 2*π
                end
                res[:weight][j] += logpdf(dist2,res[:ω_f][j]+diff)-logpdf(dist1,observed_values[i])
            elseif obs_symbol == :i
                dist1 = Cauchy(res[:i_f][j], observed_errors[i])
                dist2 = Normal(res[:i_f][j], observed_errors[i])
                res[:weight][j] += logpdf(dist2,observed_values[i])-logpdf(dist1,observed_values[i])
            elseif obs_symbol == :v_N
                dist1 = Cauchy(res[:v_N][j], observed_errors[i])
                dist2 = Normal(res[:v_N][j], observed_errors[i])
                res[:weight][j] += logpdf(dist2,observed_values[i])-logpdf(dist1,observed_values[i])
            elseif obs_symbol == :v_E
                dist1 = Cauchy(res[:v_E][j], observed_errors[i])
                dist2 = Normal(res[:v_E][j], observed_errors[i])
                res[:weight][j] += logpdf(dist2,observed_values[i])-logpdf(dist1,observed_values[i])
            elseif obs_symbol == :v_r
                dist1 = Cauchy(res[:v_r][j], observed_errors[i])
                dist2 = Normal(res[:v_r][j], observed_errors[i])
                res[:weight][j] += logpdf(dist2,observed_values[i])-logpdf(dist1,observed_values[i])
            end
        end
    end
    #adjust weights so that maximum weight is unity
    println(maximum(res[:weight]))
    res[:weight] .= res[:weight] .- maximum(res[:weight]) #weights are still in log here
    res[:weight] .= 1 #exp.(res[:weight])

    # if we did a general model, we need to weight the true anomaly
    if model_type==:general
        res[:weight] .= res[:weight].*sqrt.(1 .- res[:e_f].^2).^3 ./ (1 .+ res[:e_f].*cos.(res[:ν_i])).^2
    end

    return res
end
