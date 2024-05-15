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
function createSimpleCircularMCMCModel(observations::Vector{Symbol}, observed_values::Vector{Float64}, observed_errors::Vector{Float64};
    bhModel = arbitraryEjectaBH,
    logm1_dist::ContinuousUnivariateDistribution = Uniform(0.1,3), # in log(Msun)
    logm2_dist::ContinuousUnivariateDistribution = Uniform(0.1,3), # in log(Msun)
    logP_dist::ContinuousUnivariateDistribution = Uniform(-1,3),   # in log(days)
    vkick_dist::ContinuousUnivariateDistribution = Exponential(1), # in 100 km/s
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
        P_f = kepler_P_from_a(m1=m1, m2=m2, a=a)
        K1 = RV_semiamplitude_K1(m1=m1, m2=m2_f, P=P_f, e=e_f, i=i_f)

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
                obs_vals[i] ~ Cauchy(m1, obs_errs[i])
            elseif obs_symbol == :m2
                obs_vals[i] ~ Cauchy(m2_f, obs_errs[i])
            end
        end
    end

    return kick_model(observations, observed_values, observed_errors)

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
        θ = arccos(θ)
        xϕ ~ Normal(0,1)
        yϕ ~ Normal(0,1)
        normϕ = 1/sqrt(xϕ^2+yϕ^2)
        cosϕ = xϕ*normϕ
        ϕ = arccos(ϕ)

        #Initial systemic velocity parameters
        v_N ~ vsys_N_dist*100*km_per_s 
        v_E ~ vsys_E_dist*100*km_per_s 
        v_r ~ vsys_r_dist*100*km_per_s 

        #m1 is assumed to remain constant, no impact velocity
        a_f, e_f, Ω_f, ω_f, i_f, v_n, v_e, v_rad = 
            post_supernova_general_orbit_parameters(m1=m1, m2=m2, a=a, e=e, m2_f=m2_f, 
                vkick=vkick, θ=θ, ϕ=ϕ, ν=ν, Ω=Ω, ω=ω, i=i)
        P_f = kepler_P_from_a(m1=m1, m2=m2_f, a=a_f)
        K1 = RV_semiamplitude_K1(m1=m1, m2=m2_f, P=P_f, e=e_f, i=i_f)
        K2 = RV_semiamplitude_K1(m1=m2_f, m2=m1, P=P_f, e=e_f, i=i_f)

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
                obs_vals[i] ~ Cauchy(m1, obs_errs[i])
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
                obs_vals[i] ~ Cauchy(v_N - v_N, obs_errs[i])
            elseif obs_symbol == :v_E
                obs_vals[i] ~ Cauchy(v_E - v_E, obs_errs[i])
            elseif obs_symbol == :v_r
                obs_vals[i] ~ Cauchy(v_r - v_r, obs_errs[i])
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
    #TODO: RTW this function and the description don't match, check with Pablo

    res = Dict()
    res[:logm1] = reduce(vcat,chain[:logm1])
    res[:m1] = 10 .^ res[:logm1] *m_sun
    res[:logm2] = reduce(vcat,chain[:logm2])
    res[:m2] = 10 .^ res[:logm2] *m_sun
    res[:logP] = reduce(vcat,chain[:logP])
    res[:P] = 10 .^ res[:logP] *day
    res[:a] = kepler_a_from_P.(m1=res[:m1], m2=res[:m2], P=res[:P])
    if model_type==:simple
        res[:e] = zeros(length(res[:a]))
    elseif model_type==:general
        res[:e] = reduce(vcat,chain[:e])
        res[:ν] = reduce(vcat,shift_range.([atan(y,x) for (x,y) in zip(chain[:xν],chain[:yν])]))
        res[:Ω] = reduce(vcat,shift_range.([atan(y,x) for (x,y) in zip(chain[:xΩ],chain[:yΩ])]))
        res[:ω] = reduce(vcat,shift_range.([atan(y,x) for (x,y) in zip(chain[:xω],chain[:yω])]))
        res[:i] = reduce(vcat,[acos(x) for x in chain[:cosi]])
    else
        throw(ArgumentError("model_type=:$model_type is an invalid option"))
    end
    res[:frac] = reduce(vcat,chain[:frac])
    res[:m2_f] = bhModel.(res[:m2],res[:frac])
    res[:ejecta_mass] = res[:m2] .- res[:m2_f]
    res[:vkick] = reduce(vcat,chain[:vkick])*100*km_per_s # MCMC is done in units of 100 km/s, turn into km/s
    res[:θ] = reduce(vcat,[acos(x) for x in chain[:cosθ]])
    res[:ϕ] = reduce(vcat,shift_range.([atan(y,x) for (x,y) in zip(chain[:xϕ],chain[:yϕ])]))
    if model_type==:general
        res[:v_N] = reduce(vcat,chain[:v_N])*100*km_per_s
        res[:v_E] = reduce(vcat,chain[:v_E])*100*km_per_s
        res[:v_r] = reduce(vcat,chain[:v_r])*100*km_per_s
    end 

    sample_num = length(res[:a])
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
            post_supernova_general_orbit_parameters(m1=res[:m1][i], m2=res[:m2][i],
                a= res[:a][i], e= res[:e][i], m2_f= res[:m2_f][i], vkick= res[:vkick][i], 
                θ= res[:θ][i], ϕ= res[:ϕ][i], ν= res[:ν][i], Ω= res[:Ω][i], ω= res[:ω][i], i= res[:i][i])

            res[:a_f][i] = a_f
            res[:e_f][i] = e_f
            res[:v_N][i] = v_N + res[:v_N][i]
            res[:v_E][i] = v_E + res[:v_E][i]
            res[:v_r][i] = v_r + res[:v_r][i]
            res[:Ω_f][i] = Ω_f
            res[:ω_f][i] = ω_f
            res[:i_f][i] = i_f
        end
    elseif model_type==:simple
        res[:a_f] = Vector{Float64}(undef, sample_num)
        res[:e_f] = Vector{Float64}(undef, sample_num)
        for i in 1:sample_num
            a_f, e_f = post_supernova_circular_orbit_a(m1=res[:m1][i], m2=res[:m2][i], 
                                                         a=res[:a][i], m2_f=res[:m2_f][i], vkick=res[:vkick][i], θ=res[:θ][i], ϕ=res[:ϕ][i])
            res[:a_f][i] = a_f
            res[:e_f][i] = e_f
        end
        res[:i_f] = reduce(vcat,[acos(x) for x in chain[:cosi]])
    end
    res[:P_f] = kepler_P_from_a.(m1=res[:m1], m2=res[:m2_f], a=res[:a_f])
    res[:K1] = RV_semiamplitude_K1.(m1=res[:m1], m2=res[:m2_f], P=res[:P_f], e=res[:e_f], i=res[:i_f])
    res[:K2] = RV_semiamplitude_K1.(m1=res[:m2_f], m2=res[:m1], P=res[:P_f], e=res[:e_f], i=res[:i_f])
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
                dist1 = Cauchy(res[:m1][j], observed_errors[i])
                dist2 = Normal(res[:m1][j], observed_errors[i])
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
        res[:weight] .= res[:weight].*sqrt.(1 .- res[:e_f].^2).^3 ./ (1 .+ res[:e_f].*cos.(res[:ν])).^2
    end

    return res
end
