using LinearAlgebra

###########################################################
### 
### Standard orbital parameter calculations and conversions
### 
###########################################################

"""
    kepler_a_from_P(;m1, m2, P)

Obtain semimajor axis from period using Kepler's third law

# Arguments:
- m1: mass of first companion [g]
- m2: mass of 2nd companion   [g]
- P:  orbital period          [s]

# Output:
- a: semi-major axis of the orbit [cm]
"""
function kepler_a_from_P(;m1, m2, P)
    return @. cbrt((P)^2 *(m1+m2)*cgrav/(4.0*π^2))
end

"""
    kepler_P_from_a(;m1, m2, a)

Obtain period from semimajor axis using Kepler's third law

# Arguments:
- m1: mass of first companion      [g]
- m2: mass of 2nd companion        [g]
- a:  semi-major axis of the orbit [cm]

# Output:
- P: the orbital period            [s]
"""
function kepler_P_from_a(;m1, m2, a)
    return @. sqrt(4.0*π^2*a^3 /((m1+m2)*cgrav))
end


"""
    relative_velocity(;m1, m2, a)

Calculate the relative orbital velocity for a circular orbit.

# Arguments:
- m1: mass of first companion        [g]
- m2: mass of 2nd companion          [g]
- a:  semi-major axis of the orbit   [cm]

# Output:
- v_rel: the relative velocity [cm/s]
"""
function relative_velocity(;m1, m2, a)
    return @. sqrt(cgrav*(m1+m2)/a)
end


"""
    RV_semiamplitude_K1(;m1, m2, P, e, i)

Compute the amplitude of radial velocity variations given orbital parameters and masses

# Arguments:
- m1:   mass of observed star             [g]
- m2:   mass of companion                 [g]
- P:    orbital period                    [s]
- e:    orbital eccentricity              [-]
- i:    orbital inclination               [rad]

# Output:
- K1: amplitude of radial velocity variation of star 1  [cm/s]
"""
function RV_semiamplitude_K1(;m1, m2, P, e, i)
    return @. m2*sin(i)/sqrt(1-e^2)*cbrt(2*π*cgrav/(P*(m2 + m1)^2)) #semi amplitude of RV variation in cm/s
    
end


###########################################################
### 
### Post kick orbital properties, for circular pre-explosion orbits
### 
###########################################################

"""
    post_supernova_circular_orbit_a(;m1, m2, a, m1_f=-1, m2_f, vkick=0, θ=0, ϕ=0, vimp=0)

Compute post-kick properties for a circular pre-explosion orbit. Equivalent to
Tauris et al. (1999): Monthly Notices of the Royal Astronomical Society, Volume 310, Issue 4, pp. 1165-1169.


# Arguments:
- m1:  pre-explosion  mass of non-exploding component   [g]           
- m2:  pre-explosion  mass of exploding component       [g]       
- a:   pre-explosion orbital separation                 [cm]
- m1_f:  post-explosion mass of non-exploding component [g]           
- m2_f:  post-explosion mass of exploding component     [g]   
- vkick: kick velocity                                  [cm/s] 
- θ:     polar kick angle (away from e_par)             [rad]
- ϕ:     azimuthal kick angle (off of e_perp)           [rad]
- vimp:  imparted kick velocity on companion            [cm/s]     
    
# Output:
- a_f: post-explosion orbital separation                [cm]
- e_f: post-explosion excentricity                      [-]
"""
function post_supernova_circular_orbit_a(;m1, m2, a, m1_f=-1.0, m2_f, vkick=0.0, θ=0.0, ϕ=0.0, vimp=0.0)
    if m1_f == -1.0
        m1_f = m1
    end
    mtilde = (m1_f + m2_f)/(m1 + m2) 
    v_rel = relative_velocity(m1=m1, m2=m2, a=a)
    α = vkick/v_rel
    β = vimp/v_rel
    # convert trig functions to vars
    cosθ = cos(θ)
    sinθ = sin(θ)
    cosϕ = cos(ϕ)
    sinϕ = sin(ϕ)
    ξ = (1 + α^2 + β^2 + 2*α*cosθ - 2*α*β*sinθ*cosϕ)/(mtilde)
    if (ξ>2)
        return (NaN, NaN)
    end
    η = (α^2*sinθ^2*sinϕ^2+(1+α*cosθ)^2)/mtilde

    # Orbital parameters
    a_f = a/(2-ξ)
    e_f = sqrt(1 + (ξ-2)*η + 1e-10) # including safety floor

    return (a_f, e_f)
end


"""
    post_supernova_circular_orbit_P(;m1, m2, P, m1_f=-1, m2_f, vkick=0, θ=0, ϕ=0, vimp=0)

Same as `post_supernova_circular_orbit_a`, except that it receives the initial orbital period as
input and returns the final orbital period and eccentricity.

# Arguments:
- m1:  pre-explosion  mass of non-exploding component   [g]           
- m2:  pre-explosion  mass of exploding component       [g]       
- P:   pre-explosion orbital period                     [d]
- m1_f:  post-explosion mass of non-exploding component [g]           
- m2_f:  post-explosion mass of exploding component     [g]   
- vkick: kick velocity                                  [cm/s] 
- θ:     polar kick angle (away from e_par)             [rad]
- ϕ:     azimuthal kick angle (off of e_perp)           [rad]
- vimp:  imparted kick velocity on companion            [cm/s]     
    
# Output:
- P_f: post-explosion orbital period                    [d]
- e_f: post-explosion excentricity                      [-]
"""
function post_supernova_circular_orbit_P(;m1, m2, P, m1_f=-1.0, m2_f, vkick=0.0, θ=0.0, ϕ=0.0, vimp=0.0)
    if m1_f == -1.0
        m1_f = m1
    end
    a = kepler_a_from_P(m1=m1, m2=m2, P=P)
    (a_f, e_f) = post_supernova_circular_orbit_a(m1=m1, m2=m2, a=a,
                    m1_f=m1_f, m2_f=m2_f, vkick=vkick, θ=θ, ϕ=ϕ, vimp=vimp)
    P_f = kepler_P_from_a(m1=m1_f, m2=m2_f, a=a_f)

    return (P_f,e_f)
end

"""
    post_supernova_circular_orbit_vsys(;m1, m2, a, m1_f=-1, m2_f, vkick=0, θ=0, ϕ=0, vimp=0)

Compute post-kick properties for a circular pre-explosion orbit using equations from Tauris & Takens (1999)

# Arguments:
- m1:  pre-explosion  mass of non-exploding component   [g]           
- m2:  pre-explosion  mass of exploding component       [g]       
- a:   pre-explosion orbital separation                 [cm]
- m1_f:  post-explosion mass of non-exploding component [g]           
- m2_f:  post-explosion mass of exploding component     [g]   
- vkick: kick velocity                                  [cm/s] 
- θ:     polar kick angle (away from e_par)             [rad]
- ϕ:     azimuthal kick angle (off of e_perp)           [rad]
- vimp:  imparted kick velocity on companion            [cm/s]     
    
# Output:
- vsys_f: post-explosion systemic velocity              [cm/s]
"""
function post_supernova_circular_orbit_vsys(;m1, m2, a, m1_f=-1, m2_f, vkick=0, θ=0, ϕ=0, vimp=0)
    if m1_f == -1
        m1_f = m1
    end
    v_rel = relative_velocity(m1=m1, m2=m2, a=a)
    # convert trig functions to vars
    cosθ = cos(θ)
    sinθ = sin(θ)
    cosϕ = cos(ϕ)
    sinϕ = sin(ϕ)

    # Systemic velocity
    Δp_x = (m2_f*m1 - m2*m1_f)/(m2 + m1)*v_rel + m2_f*vkick*cosθ
    Δp_y = m1_f*vimp + m2_f*vkick*sinθ*cosϕ
    Δp_z = m2_f*vkick*sinθ*sinϕ
    vsys_f = sqrt(Δp_x^2 + Δp_y^2 + Δp_z^2)/(m2_f + m1_f)

    return vsys_f
end



###########################################################
### 
### Post kick orbital properties, for general pre-explosion orbits
### 
###########################################################

"""
    post_supernova_general_orbit_parameters(;m1, m2, a, e=0, m1_f=-1, m2_f, vkick=0, θ=0, ϕ=0, vimp=0,
        ν=0, Ω=0, ω=0, i=0)

Compute post-kick properties for a general pre-explosion orbit 
using equations from [Marchant, Willcox, Vigna-Gomez] TODO

# Arguments:
- m1:  pre-explosion  mass of non-exploding component    [g]           
- m2:  pre-explosion  mass of exploding component        [g]       
- a:   pre-explosion orbital separation                  [cm]
- e:   pre-explosion orbital eccentricity                [-]
- m1_f:  post-explosion mass of non-exploding component  [g]           
- m2_f:  post-explosion mass of exploding component      [g]   
-
- vkick: kick velocity                                   [cm/s] 
- θ:     polar kick angle (away from e_par)              [rad]
- ϕ:     azimuthal kick angle (off of e_perp)            [rad]
- vimp:  imparted kick velocity on companion             [cm/s]     
- Initial orbital orientation angles: 
    - ν: true anomaly                                    [rad]
    - Ω: pre-explosion longitude of the ascending node   [rad]
    - ω: pre-explosion argument of periastron            [rad]
    - i: pre-explosion inclination                       [rad]
    
# Output: RTW: check! 
- a_f:   post-explosion orbital separation               [cm]
- e_f:   post-explosion orbital eccentricity             [-]
- Ω_f:   post-explosion longitude of ascending node      [rad]      
- ω_f:   post-explosion argument of periastron           [rad]    
- i_f:   post-explosion inclination                      [rad]     
- v_n:   post-explosion systemic velocity, toward N      [rad]
- v_w:   post-explosion systemic velocity, toward W      [rad]      
- v_rad: post-explosion radial velocity, toward O        [rad]      
"""
function post_supernova_general_orbit_parameters(;m1, m2, a, e=0, m1_f=-1, m2_f, vkick=0, 
        θ=0, ϕ=0, vimp=0, ν=0, Ω=0, ω=0, i=0)
    if m1_f == -1
        m1_f = m1
    end
    M = m1 + m2
    M_f = m1_f + m2_f

    # convert trig functions to vars
    cosθ = cos(θ)
    sinθ = sin(θ)
    cosϕ = cos(ϕ)
    sinϕ = sin(ϕ)
    cosν = cos(ν)
    sinν = sin(ν)
    cosi = cos(i)
    sini = sin(i)
    cosω = cos(ω)
    sinω = sin(ω)
    cosΩ = cos(Ω)
    sinΩ = sin(Ω)

    # construct useful intermediary parameters
    f_ν = (1 - e^2)/(1 + e*cosν)
    g_ν = sqrt((1 + 2*e*cosν + e^2)/(1 - e^2))
    h_ν = -e*sinν/sqrt(1 + 2*e*cosν + e^2)
    j_ν = (1 + e*cosν)/sqrt(1 + 2*e*cosν + e^2)

    v_rel = g_ν*sqrt(cgrav*M/a)
    α = vkick/v_rel
    β = vimp/v_rel
    ξ = f_ν*g_ν^2*M/M_f* (1 + α^2 + β^2 + 2* (α*cosθ - h_ν*β* (1 + α*cosθ) - j_ν*β*α*sinθ*cosϕ))
    if ξ>2
        return (NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
    end

    a_f = f_ν*a/(2 - ξ)
    Lvec_norm = sqrt(α^2*sinθ^2*sinϕ^2 + (h_ν*α*sinθ*cosϕ - j_ν*(1 + α*cosθ))^2)
    η = f_ν*g_ν^2*M/M_f*Lvec_norm^2
    e_f = sqrt(1 + (ξ-2)*η + 1e-10) # including safety floor

    #angle between x and direction of motion
    # perhaps can be skipped to save some time, just get directly cos and sin
    τ = acos(max(-1, min(1, j_ν)))
    if (sinν < 0)
        τ = -τ
    end
    #δ = ω + ν - τ
    cosτ = cos(τ)
    sinτ = sin(τ)
    cosδ = cosν*sinτ*sinω - sinν*cosτ*sinω + sinν*sinτ*cosω + cosν*cosτ*cosω
    sinδ = sinν*sinτ*sinω + cosν*cosτ*sinω - cosν*sinτ*cosω + sinν*cosτ*cosω

    # Elements of rotation matrix
    R_w_par =  cosi*cosΩ*cosδ - sinΩ*sinδ
    R_w_per = -cosi*cosΩ*sinδ - sinΩ*cosδ
    R_w_z   = -sini*cosΩ
    R_n_par =  cosi*sinΩ*cosδ + cosΩ*sinδ 
    R_n_per = -cosi*sinΩ*sinδ + cosΩ*cosδ 
    R_n_z   = -sini*sinΩ
    R_o_par =  sini*cosδ
    R_o_per = -sini*sinδ
    R_o_z   =  cosi

    # velocity, simply compute from change in momentum
    v_par = (-(m2 - m2_f)*m1/M*v_rel + (m1 - m1_f)*m2/M*v_rel
             + m2_f*vkick*cosθ + h_ν*m1_f*vimp)/M_f
    v_per = (m2_f*vkick*sinθ*cosϕ + j_ν*m1_f*vimp)/M_f
    v_z   = m2_f*vkick*sinθ*sinϕ/M_f

    v_w = R_w_par*v_par + R_w_per*v_per + R_w_z*v_z
    v_n = R_n_par*v_par + R_n_per*v_per + R_n_z*v_z
    v_o = R_o_par*v_par + R_o_per*v_per + R_o_z*v_z
    v_rad = -v_o # swap because v_o is defined towards the observer, but we're interested in the opposite

    # obtain inclination from direction of orbital angular momentum vector
    L_par = j_ν*α*sinθ*sinϕ/Lvec_norm
    L_per = -h_ν*α*sinθ*sinϕ/Lvec_norm
    L_z = h_ν*α*sinθ*cosϕ - j_ν*(1 + α*cosθ)/Lvec_norm

    L_w = R_w_par*L_par + R_w_per*L_per + R_w_z*L_z
    L_n = R_n_par*L_par + R_n_per*L_per + R_n_z*L_z
    L_o = R_o_par*L_par + R_o_per*L_per + R_o_z*L_z

    i_f = acos(min(1, abs(L_o)))

    # compute longitude of the ascending node
    # this is obtained by taking the vector n = k x L, where k=(0,0,1) points to the observer
    # then n = (n_w, n_n, 0) points to the ascending node
    # TODO RTW: need to check the code below, the subscripts look wrong
    n_w = -L_n
    n_n = L_w
    n_norm = sqrt(n_w^2 + n_n^2)
    Ω_f = acos(max(-1, min(1, n_n/n_norm)))
    if n_w>0
        Ω_f = 2π - Ω_f
    end
    # compute post-explosion argument of periastron
    if (e_f > 0)
        periastron_angle = acos(max(-1, min(1, 1/e_f*(a_f/(a*f_ν)*(1 - e_f^2)-1))))
    else
        periastron_angle = 0 # need to check this, though maybe irrelevant as chance of this is null
    end
    # the periastron angle is the same as the true anomaly if the star is moving
    # away from periastron. We need to compute the component of velocity
    # along the line joining both objects (in the COM frame).
    v1y_cm = vimp + h_ν*(m2/M*v_rel - v_par) - j_ν*v_per
    if v1y_cm>0
        periastron_angle = 2π-periastron_angle
    end
    # compute the angle from current location to ascending node
    # in the rotated frame. This is R*\hat{y}
    rvec_w = R_w_par*h_ν + R_w_per*j_ν
    rvec_n = R_n_par*h_ν + R_n_per*j_ν
    rvec_o = R_o_par*h_ν + R_o_per*j_ν
    dot_prod = (n_w*rvec_w + n_n*rvec_n)/n_norm
    angle_to_asc_node = acos(max(-1, min(1, dot_prod)))
    # We need to know if ν is before or after the ascending node.
    # To do this we take the cross product of the vector pointing from the COM
    # to star 1, and that pointing to the ascending node. If this points in
    # the direction of the orbital angular momentum, then we have not
    # overtaken the ascending node
    cross_vec_w = -rvec_o*n_n
    cross_vec_n = rvec_o*n_w
    cross_vec_o = rvec_w*n_n - rvec_n*n_w
    cross_vec_dot_L = cross_vec_w*L_w + cross_vec_n*L_n + cross_vec_o*L_o

    if cross_vec_dot_L > 0
        ω_f = angle_to_asc_node - periastron_angle
    else
        ω_f = 2π - (angle_to_asc_node + periastron_angle)
    end
    if ω_f < 0
        ω_f = 2π + ω_f
    end

    return (a_f, e_f, Ω_f, ω_f, i_f, v_n, v_w, v_rad)
end




