using LinearAlgebra
#using Symbolics
#using IfElse

###########################################################
### 
### Standard orbital parameter calculations and conversions
### 
###########################################################

"""
    kepler_a_from_P(P, m1, m2)

Obtain semimajor axis from period using Kepler's third law

# Arguments:
- P:  orbital period          [d]
- m1: mass of first companion [Msun]
- m2: mass of 2nd companion   [Msun]

# Output:
- a: semi-major axis of the orbit [cm]
"""
function kepler_a_from_P(P, m1, m2)
    return cbrt((P*day_in_sec)^2*cgrav*(m1+m2)*m_sun/(4.0*π^2))
end


"""
    kepler_P_from_a(a, m1, m2)

Obtain period from semimajor axis using Kepler's third law

# Arguments:
- a:  semi-major axis of the orbit [cm]
- m1: mass of first companion      [Msun]
- m2: mass of 2nd companion        [Msun]

# Output:
- P: the orbital period            [d]
"""
function kepler_P_from_a(a, m1, m2)
    return sqrt(4.0*π^2*a^3/(cgrav*(m1+m2)*m_sun))/(day_in_sec)
end

"""
    RV_semiamplitude_K1(P, e, i, m1, m2)

Compute the amplitude of radial velocity variations given orbital parameters and masses

# Arguments:
- P:     orbital period                    [d]
- e:     orbital eccentricity              [dimensionless]
- i:     orbital inclination               [rad]
- m1:    mass of observed star             [Msun]
- m2:    mass of companion                 [Msun]

# Output:
- K1: amplitude of radial velocity variation of star 1  [km/s]
"""
function RV_semiamplitude_K1(P, e, i, m1, m2)
    # RTW where does the 1e5 come from here?
    return (m2*sin(i))*cbrt(2*π*cgrav/(P*day_in_sec)*m_sun/(m2+m1)^2)/sqrt(1-e^2)/1e5 #semi amplitude of RV variation in km/s
end


###########################################################
### 
### Post kick orbital properties, for circular pre-explosion orbits
### 
###########################################################

"""
    post_circular_kick_P_e_a(m1i, m2i, m2f, vrel, vkick, theta, phi; [Pi, ai])

Compute post-kick properties for a circular pre-explosion orbit using equations from Tauris et al. (1999)

# Arguments:
- m1i:   mass of non-exploding component             [Msun]           
- m2i:   mass of exploding component                 [Msun]       
- m2f:   post explosion mass of exploding component  [Msun]   
- vrel:  pre-explosion relative velocity             [km/s] 
- vkick: kick velocity                               [km/s] 
- θ: polar kick angle (away from e_par)              [rad]
- ϕ: azimuthal kick angle (off of e_perp)            [rad]
- Only one of:
    - Pi:    pre-explosion orbital period            [d]
    - ai:    pre-explosion orbital separation        [cm]
    
# Output:
- Pf: post-explosion orbital period                  [d]
- ef: post-explosion excentricity                    [dimensionless]
- af: post-explosion orbital separation              [cm]
- vsysf: post-explosion systemic velocity            [km/s]

## RTW why does m1i stay fixed here, but is allowed to change in the non-circular case?
"""
function post_circular_kick_parameters(m1i, m2i, m2f, vrel, vkick, θ, ϕ; Pi=nothing, ai=nothing)
    # Exactly one of P or a is required
    if (isnothing(Pi) == isnothing(ai))
        throw(ArgumentError("Exactly one of Pi or ai is required"))
    elseif isnothing(Pi)
        Pi = kepler_P_from_a(ai, m1i, m2i)
    elseif isnothing(ai)
        ai = kepler_a_from_P(Pi, m1i, m2i)
    end
    mtilde = (m1i+m2f)/(m1i+m2i) 
    vkick_div_vrel = vkick/vrel
    # convert trig functions to vars
    cosθ = cos(θ)
    sinθ = sin(θ)
    cosϕ = cos(ϕ)
    sinϕ = sin(ϕ)
    xi = (1+vkick_div_vrel^2+2*vkick_div_vrel*cosθ)/(mtilde)
    Q = xi - 1 - (vkick_div_vrel*sinθ*cosϕ)^2/(mtilde)

    if (xi>2)
        return (NaN, NaN)
    end

    # Orbital parameters
    Pf = Pi/sqrt(mtilde*(2-xi)^3)
    ef = sqrt(1+(xi-2)*(Q+1))
    af = ai/(2-xi)

    # Systemic velocity
    Δp_x = (m2f*m1i-m2i*m1i)/(m2i+m1i)*vrel+m2f*vkick*cosθ
    Δp_y = m2f*vkick*sinθ*cosϕ
    Δp_z = m2f*vkick*sinθ*sinϕ
    vsysf = sqrt(Δp_x^2+Δp_y^2+Δp_z^2)/(m2f+m1i)/1e5 
    # RTW - where does 1e5 come from?

    return (Pf, ef, af, vsysf)
end

###########################################################
### 
### Post kick orbital properties, for general pre-explosion orbits
### 
###########################################################

"""
    post_general_kick_parameters(m1i, m2i, m1f, m2f, ei, vkick, vimp, 
    ν, θ, ϕ, Ω, ω, ι; [Pi, ai])

Compute post-kick properties for a circular pre-explosion orbit using equations from Tauris et al. (1999)

# Arguments:
- m1i:   mass of non-exploding component             [Msun]           
- m2i:   mass of exploding component                 [Msun]       
- m2f:   post-explosion mass of exploding component  [Msun]   
- ei:    pre-explosion orbital eccentricity          [dimensionless]
- vkick: kick velocity on exploding component        [km/s] 
- vimp:  imparted kick velocity on companion         [km/s]     
- Kick and orientation angles: 
    - θ: polar kick angle (away from e_par)          [rad]
    - ϕ: azimuthal kick angle (off of e_perp)        [rad]
    - νi: true anomaly                                [rad]
    - Ωi: pre-explosion longitude of the ascending node
    - ωi: pre-explosion argument of periastron
    - ιi: pre-explosion inclination
- Only one of:
    - Pi: pre-explosion orbital period               [d]
    - ai: pre-explosion orbital separation           [cm]
    
# Output: RTW: check! 
- af:   post-explosion orbital separation            [cm]
- ef:   post-explosion orbital eccentricity          [dimensionless]
- Ωf:   post-explosion longitude of ascending node   [rad]      
- ωf:   post-explosion argument of periastron        [rad]    
- ιf:   post-explosion inclination                   [rad]     
- v_n:  post-explosion systemic velocity, toward N   [rad]
- v_e:  post-explosion systemic velocity, toward E   [rad]      
- -v_o: post-explosion systemic velocity, toward O   [rad]      
"""
function post_general_kick_parameters(m1i, m2i, m1f, m2f, ei, vkick, vimp, 
    θ, ϕ, νi, Ωi, ωi, ιi; Pi=nothing, ai=nothing)

    # Exactly one of P or a is required
    if (isnothing(Pi) == isnothing(ai))
        throw(ArgumentError("Exactly one of Pi or ai is required"))
    elseif isnothing(Pi)
        Pi = kepler_P_from_a(ai, m1i, m2i)
    elseif isnothing(ai)
        ai = kepler_a_from_P(Pi, m1i, m2i)
    end

    # convert trig functions to vars
    cosθ = cos(θ)
    sinθ = sin(θ)
    cosϕ = cos(ϕ)
    sinϕ = sin(ϕ)
    cosνi = cos(νi)
    sinνi = sin(νi)
    cosιi = cos(ιi)
    sinιi = sin(ιi)
    cosωi = cos(ωi)
    sinωi = sin(ωi)
    cosΩi = cos(Ωi)
    sinΩi = sin(Ωi)

    # construct useful intermediary parameters
    f_ν = (1-ei^2)/(1+ei*cosνi)
    g_ν = sqrt((1+2*ei*cosνi+ei^2)/(1-ei^2))
    Mi = m1i+m2i
    Mf = m1f+m2f
    vrel = g_ν*sqrt(cgrav*Mi/ai)

    h_ν = -ei*sinνi/sqrt(1+2*ei*cosνi+ei^2)
    j_ν = (1+ei*cosνi)/sqrt(1+2*ei*cosνi+ei^2)

    α = vkick/vrel
    β = vimp/vrel

    ξ = f_ν*g_ν^2*Mi/Mf*(1+α^2+β^2+2*(-h_ν*β*(1+α*cosθ)+α*cosθ-j_ν*β*α*sinθ*cosϕ))

    if ξ>2
        return (NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
    end

    af = f_ν*ai/(2-ξ)

    Lvec_norm = sqrt(α^2*sinθ^2*sinϕ^2+(h_ν*α*sinθ*cosϕ-j_ν*(1+α*cosθ))^2)
    η = f_ν*g_ν^2*Mi/Mf*Lvec_norm^2

    ef = sqrt(1+(ξ-2)*η)

    #angle between x and direction of motion
    # perhaps can be skipped to save some time, just get directly cos and sin
    τ = acos(max(-1, min(1, j_ν)))
    if (sinνi < 0)
        τ = -τ
    end
    #δ = ω + ν - τ
    cosτ = cos(τ)
    sinτ = sin(τ)
    cosδ = cosνi*sinτ*sinωi - sinνi*cosτ*sinωi + sinνi*sinτ*cosωi + cosνi*cosτ*cosωi
    sinδ = sinνi*sinτ*sinωi + cosνi*cosτ*sinωi - cosνi*sinτ*cosωi + sinνi*cosτ*cosωi

    # Elements of rotation matrix
    # RTW: this coordinate system uses E, N, and O, but the text refers instead to celestial W, not E... need to double check
    R_e_par = cosιi*cosΩi*cosδ-sinΩi*sinδ
    R_e_per = -(cosιi*cosΩi*sinδ)-sinΩi*cosδ
    R_e_z   = -(sinιi*cosΩi)
    R_n_par = cosΩi*sinδ+cosιi*sinΩi*cosδ
    R_n_per = cosΩi*cosδ-cosιi*sinΩi*sinδ
    R_n_z   = -(sinιi*sinΩi)
    R_o_par = sinιi*cosδ
    R_o_per = -(sinιi*sinδ)
    R_o_z   = cosιi

    # velocity, simply compute from change in momentum
    v_par = (-(m2i-m2f)*m1i/Mi*vrel +(m1i-m1f)*m2i/Mi*vrel
                    + m2f*vkick*cosθ + h_ν*m1f*vimp)/Mf
    v_per = (m2f*vkick*sinθ*cosϕ+j_ν*m1f*vimp)/Mf
    v_z = (m2f*vkick*sinθ*sinϕ)/Mf

    v_e = R_e_par*v_par + R_e_per*v_per + R_e_z*v_z
    v_n = R_n_par*v_par + R_n_per*v_per + R_n_z*v_z
    v_o = R_o_par*v_par + R_o_per*v_per + R_o_z*v_z

    # obtain inclination from direction of orbital angular momentum vector
    L_par = j_ν*α*sinθ*sinϕ/Lvec_norm
    L_per = -h_ν*α*sinθ*sinϕ/Lvec_norm
    L_z = h_ν*α*sinθ*cosϕ-j_ν*(1+α*cosθ)/Lvec_norm

    L_e = R_e_par*L_par + R_e_per*L_per + R_e_z*L_z
    L_n = R_n_par*L_par + R_n_per*L_per + R_n_z*L_z
    L_o = R_o_par*L_par + R_o_per*L_per + R_o_z*L_z

    ιf = acos(min(1, abs(L_o)))

    # compute longitude of the ascending node
    # this is obtained by taking the vector n = k x L, where k=(0,0,1) points to the observer
    # then n = (n_e, n_n, 0) points to the ascending node
    n_e = -L_n
    n_n = L_e
    n_norm = sqrt(n_e^2+n_n^2)
    Ωf = acos(max(-1, min(1, n_n/n_norm)))
    if n_e>0
        Ωf = 2π - Ωf
    end
    # compute post-explosion argument of periastron
    if (ef > 0)
        periastron_angle = acos(max(-1, min(1, 1/ef*(af/(a*f_ν)*(1-ef^2)-1))))
    else
        periastron_angle = 0 # need to check this, though maybe irrelevant as chance of this is null
    end
    # the periastron angle is the same as the true anomaly if the star is moving
    # away from periastron. We need to compute the component of velocity
    # amont the line joining both objects (in the COM frame).
    v1y_cm = vimp + h_ν*(m2/M*vrel - v_par) - j_ν*v_per
    if v1y_cm>0
        periastron_angle = 2π-periastron_angle
    end
    # compute the angle from current location to ascending node
    # in the rotated frame. This is R*\hat{y}
    rvec_e = R_e_par*h_ν + R_e_per*j_ν
    rvec_n = R_n_par*h_ν + R_n_per*j_ν
    rvec_o = R_o_par*h_ν + R_o_per*j_ν
    dot_prod = (n_e*rvec_e + n_n*rvec_n)/n_norm
    angle_to_asc_node = acos(max(-1, min(1, dot_prod)))
    # We need to know if ν is before or after the ascending node.
    # To do this we take the cross product of the vector pointing from the COM
    # to star 1, and that pointing to the ascending node. If this points in
    # the direction of the orbital angular momentum, then we have not
    # overtaken the ascending node
    cross_vec_e = - rvec_o*n_n
    cross_vec_n = rvec_o*n_e
    cross_vec_o = rvec_e*n_n - rvec_n*n_e
    cross_vec_dot_L = cross_vec_e*L_e + cross_vec_n*L_n + cross_vec_o*L_o

    if cross_vec_dot_L > 0
        ωf = angle_to_asc_node-periastron_angle
    else
        ωf = 2π - (angle_to_asc_node+periastron_angle)
    end
    if ωf < 0
        ωf = 2π + ωf
    end

    # RTW why do we return -v_o instead of positive?
    return (af, ef, Ωf, ωf, ιf, v_n, v_e, -v_o)
end

