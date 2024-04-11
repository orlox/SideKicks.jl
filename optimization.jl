using Turing
using StatsPlots
# solar mass and gravitational constant in CGS
const cgrav = 6.67430e-8
const mu_sun = 1.3271244e26
const m_sun = mu_sun/cgrav
const day_in_sec = 24*3600

##
"""
    post_kick_parameters_P(P,m1,m2i,m2f,vkick)

Compute post-kick properties using equations from Tauris et al. (1999), but using the orbital period

#Arguments:
- P:     pre-explosion orbital period               [d]
- mtilde: (m1_i+m2_f)/(m1_i+m2_i) [unitless], here we have
    - m1:    mass of non-exploding component
    - m2i:   mass of exploding component
    - m2f:   post explosion mass of exploding component
- vkick_div_vrel: kick velocity in units of the pre-explosion relative velocity    [unitless]
    - vrel = sqrt(G(m1_i + m2_i)/a)
- theta: polar angle indicating kick direction      [rad]
- phi:   azimuthal angle indicating kick direction  [rad]

#Output: (Ppost, epost) where
- Ppost: post-explosion orbital period              [d]
- epost: post-explosion eccentricity                [dimensionless]
"""
function post_kick_parameters_P(P,mtilde,vkick_div_vrel,theta,phi)
    xi = (1+vkick_div_vrel^2+2*vkick_div_vrel*cos(theta))/(mtilde)
    Q = xi - 1 - (vkick_div_vrel*sin(theta)*cos(phi))^2/(mtilde)

    if (xi>2)
        return (NaN*P, NaN*P) # multiply by P to make outcome type stable
    end
    Pfinal = P/sqrt(mtilde*(2-xi)^3)
    efinal = sqrt(1+(xi-2)*(Q+1))
    return (Pfinal, efinal)
end

##

typeof(post_kick_parameters_P)
first(methodswith(Function), 5)
##

"""
    RV_semiamplitude_K(P, e, i, m1, m2)

Compute the amplitude of radial velocity variations given orbital parameters and masses

#Arguments:
- P:     orbital period                    [d]
- e:     orbital eccentricity              [dimensionless]
- i:     orbital inclination               [rad]
- m1:    mass of observed star             [Msun]
- m2:    mass of companion                 [Msun]

#Output:
- K1: amplitude of radial velocity variation of star 1  [km/s]
"""
function RV_semiamplitude_K(P, e, i, m1, m2)
    return (m2*sin(i))*cbrt(2*π*cgrav/(P*day_in_sec)*m_sun/(m2+m1)^2)/sqrt(1-e^2)/1e5 #semi amplitude of RV variation in km/s
end

##
# try to find models that mathc a given period and eccentricity post explosion
@model function kick_model(Pobs, Pobs_err, eobs, eobs_err, M1obs, M1obs_err, K1obs, K1obs_err)
    logm1_i ~ Uniform(0,2)
    logm2_i ~ Uniform(0,2)
    logP_i ~ Uniform(0,2)
    vkick ~ Exponential(1) #in 100 km/s
    frac_lost ~ Uniform(0,1) #fraction of mass lost from exploding star
    cosθ ~ Uniform(-1,1)
    θ = acos(cosθ)
    xϕ ~ Normal(0,1)
    yϕ ~ Normal(0,1)
    cosi ~ Uniform(0,1) # random inclination for post-explosion orbit
    i_f = acos(cosi)
    normϕ = 1/sqrt(xϕ^2+yϕ^2)
    cosϕ = xϕ*normϕ
    sinϕ = yϕ*normϕ
    ϕ = acos(cosϕ)
    if sinϕ < 0
        ϕ = 2π - ϕ
    end

    m1_i = 10^logm1_i
    m2_i = 10^logm2_i
    m2_f = (1-frac_lost)*m2_i
    P_i = 10^logP_i
    #get separation from Kepler's third law
    a_i = cbrt((P_i*day_in_sec)^2*cgrav*(m1_i+m2_i)*m_sun/(4.0*π^2))
    vrel = sqrt(cgrav*(m1_i+m2_i)*m_sun/a_i)
    vkick_div_vrel = vkick*1e7/vrel
    mtilde = (m1_i+m2_f)/(m1_i+m2_i)

    Pfinal, efinal = post_kick_parameters_P(P_i,mtilde,vkick_div_vrel,θ,ϕ)

    K1 = RV_semiamplitude_K(Pfinal, efinal, i_f, m1_i, m2_f)

    Pobs ~ Normal(Pfinal, Pobs_err)
    eobs ~ Normal(efinal, eobs_err)
    M1obs ~ Normal(m1_i, M1obs_err)
    K1obs ~ Normal(K1, K1obs_err)
end

##
model = kick_model(10.4031, 0.01, 0.017, 0.012, 25.0, 2.3, 81.4, 1.3)
iterations = 10_000
@time chain = sample(model, NUTS(10_000,0.8), iterations);

##
plot(chain)