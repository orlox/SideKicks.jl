using LinearAlgebra

###########################################################
### 
### Standard orbital parameter calculations and conversions
### 
###########################################################

"""
    kepler_a_from_P(;m_1, m_2, P)

Obtain semimajor axis from period using Kepler's third law

# Arguments:
- m_1: mass of first companion [Msun]
- m_2: mass of 2nd companion   [Msun]
- P:  orbital period           [d]

# Output:
- a: semi-major axis of the orbit [cm]
"""
function kepler_a_from_P(;m_1, m_2, P)
    return cbrt((P*day_in_sec)^2*cgrav*(m_1+m_2)*m_sun/(4.0*π^2))
end

"""
    kepler_P_from_a(;m_1, m_2, a)

Obtain period from semimajor axis using Kepler's third law

# Arguments:
- m_1: mass of first companion      [Msun]
- m_2: mass of 2nd companion        [Msun]
- a:  semi-major axis of the orbit  [cm]

# Output:
- P: the orbital period            [d]
"""
function kepler_P_from_a(;m_1, m_2, a)
    return sqrt(4.0*π^2*a^3/(cgrav*(m_1+m_2)*m_sun))/(day_in_sec)
end


"""
    average_orbital_velocity(;m_1, m_2, a)

Calculate the average orbital velocity, or equivalently
the fixed orbital velocity for a circular orbit.

# Arguments:
- m_1: mass of first companion        [Msun]
- m_2: mass of 2nd companion          [Msun]
- a:  semi-major axis of the orbit    [cm]

# Output:
- v_orb: the average orbital velocity [km/s]
"""
function average_orbital_velocity(;m_1, m_2, a)
    return sqrt(cgrav*(m_1+m_2)*m_sun/ a)/(km) 
end


"""
    RV_semiamplitude_K1(;m_1, m_2, P, e, i)

Compute the amplitude of radial velocity variations given orbital parameters and masses

# Arguments:
- m_1:   mass of observed star             [Msun]
- m_2:   mass of companion                 [Msun]
- P:     orbital period                    [d]
- e:     orbital eccentricity              [dimensionless]
- i:     orbital inclination               [rad]

# Output:
- K1: amplitude of radial velocity variation of star 1  [km/s]
"""
function RV_semiamplitude_K1(;m_1, m_2, P, e, i)
    return (m_2*sin(i))*cbrt(2*π*cgrav/(P*day_in_sec)*m_sun/(m_2+m_1)^2)/sqrt(1-e^2)/(km) #semi amplitude of RV variation in km/s
end


###########################################################
### 
### Post kick orbital properties, for circular pre-explosion orbits
### 
###########################################################

"""
    post_supernova_circular_orbit_a(;m_1i, m_2i, a_i, m_1f=-1, m_2f, vkick=0, θ=0, ϕ=0, vimp=0)

Compute post-kick properties for a circular pre-explosion orbit using equations 
from Tauris & Takens (1999), except that here we have M2 as the exploding star.

# Arguments:
- m_1i:  pre-explosion  mass of non-exploding component [Msun]           
- m_2i:  pre-explosion  mass of exploding component     [Msun]       
- a_i:   pre-explosion orbital separation               [cm]
- m_1f:  post-explosion mass of non-exploding component [Msun]           
- m_2f:  post-explosion mass of exploding component     [Msun]   
- vkick: kick velocity                                  [km/s] 
- θ:     polar kick angle (away from e_par)             [rad]
- ϕ:     azimuthal kick angle (off of e_perp)           [rad]
- vimp:  imparted kick velocity on companion            [km/s]     
    
# Output:
- a_f: post-explosion orbital separation                [cm]
- e_f: post-explosion excentricity                      [dimensionless]
"""
function post_supernova_circular_orbit_a(;m_1i, m_2i, a_i, m_1f=-1, m_2f, vkick=0, θ=0, ϕ=0, vimp=0)
    if m_1f == -1
        m_1f = m_1i
    end
    mtilde = (m_1f+m_2f)/(m_1i+m_2i) 
    v_rel = average_orbital_velocity(;m_1=m_1i, m_2=m_2i, a=a_i)
    α = vkick/v_rel
    β = vimp/v_rel
    #vkick_div_vrel = vkick/v_rel
    #vimp_div_vrel = vimp/v_rel
    # convert trig functions to vars
    cosθ = cos(θ)
    sinθ = sin(θ)
    cosϕ = cos(ϕ)
    sinϕ = sin(ϕ)

    #println("==> ", cosθ, " ", cosϕ)
    ξ = (1 + α^2 + β^2 + 2*α*cosθ - 2*β*sinθ*cosϕ)/(mtilde)
    #println(ξ)
    #ξ = f_ν*g_ν^2*M_i/M_f*(1+α^2+β^2+2*(-h_ν*β*(1+α*cosθ)+α*cosθ-j_ν*β*α*sinθ*cosϕ))



    if (ξ>2)
        return (NaN, NaN)
    end
    tanγ = (vkick*sinθ*cosϕ - vimp)/sqrt((v_rel + vkick*cosθ)^2 + vkick^2*sinθ^2*sinϕ^2)
    cosγ = 1/sqrt(1+tanγ^2)

    # Orbital parameters
    a_f = a_i/(2-ξ)
    e_f = sqrt(1 + ξ*(ξ-2)*cosγ^2)
    #println("af=",a_f)

    #Q = (ξ-1) - (vkick_div_vrel*sinθ*cosϕ - vimp_div_vrel)^2/(mtilde) 
    #e_f_old = sqrt(1+(ξ-2)*(Q+1)) # RTW where does this come from?
    # RTW TODO test both values for e_f

    return (a_f, e_f)
end


"""
    post_supernova_circular_orbit_P(;m_1i, m_2i, P_i, m_1f=-1, m_2f, vkick=0, θ=0, ϕ=0, vimp=0)

Compute post-kick properties for a circular pre-explosion orbit using equations 
from Tauris & Takens (1999), except that here we have M2 as the exploding star.

Note: This function should be slightly slower than post_supernova_circular_orbit_a, thus
that one is preferred over this one.

# Arguments:
- m_1i:  pre-explosion  mass of non-exploding component [Msun]           
- m_2i:  pre-explosion  mass of exploding component     [Msun]       
- P_i:   pre-explosion orbital period                   [d]
- m_1f:  post-explosion mass of non-exploding component [Msun]           
- m_2f:  post-explosion mass of exploding component     [Msun]   
- vkick: kick velocity                                  [km/s] 
- θ:     polar kick angle (away from e_par)             [rad]
- ϕ:     azimuthal kick angle (off of e_perp)           [rad]
- vimp:  imparted kick velocity on companion            [km/s]     
    
# Output:
- P_f: post-explosion orbital period                    [d]
- e_f: post-explosion excentricity                      [dimensionless]
"""
function post_supernova_circular_orbit_P(;m_1i, m_2i, P_i, m_1f=-1, m_2f, vkick=0, θ=0, ϕ=0, vimp=0)
    if m_1f == -1
        m_1f = m_1i
    end
    mtilde = (m_1f+m_2f)/(m_1i+m_2i) 
    a_i = kepler_a_from_P(;m_1=m_1i, m_2=m_2i, P=P_i)
    v_rel = average_orbital_velocity(;m_1=m_1i, m_2=m_2i, a=a_i)
    vkick_div_vrel = vkick/v_rel
    vimp_div_vrel = vimp/v_rel
    # convert trig functions to vars
    cosθ = cos(θ)
    sinθ = sin(θ)
    cosϕ = cos(ϕ)
    sinϕ = sin(ϕ)
    ξ = (1 + vkick_div_vrel^2 + vimp_div_vrel^2 + 2*vkick_div_vrel*cosθ - 2*vimp_div_vrel*sinθ*cosϕ)/(mtilde)
    if (ξ>2)
        return (NaN, NaN)
    end
    tanγ = (vkick*sinθ*cosϕ - vimp)/sqrt((v_rel + vkick*cosθ)^2 + vkick^2*sinθ^2*sinϕ^2)
    cosγ = 1/sqrt(1+tanγ^2)

    # Orbital parameters
    P_f = P_i/sqrt(mtilde*(2-ξ)^3) # RTW check this
    e_f = sqrt(1 + ξ*(ξ-2)*cosγ^2)

    #Q = (ξ-1) - (vkick_div_vrel*sinθ*cosϕ - vimp_div_vrel)^2/(mtilde) 
    #e_f_old = sqrt(1+(ξ-2)*(Q+1)) # RTW where does this come from?
    # RTW TODO test both values for e_f

    return (P_f, e_f)
end


"""
    post_supernova_circular_orbit_vsys(;m_1i, m_2i, a_i, m_1f=-1, m_2f, vkick=0, θ=0, ϕ=0, vimp=0)

Compute post-kick properties for a circular pre-explosion orbit using equations from Tauris & Takens (1999)

# Arguments:
- m_1i:  pre-explosion  mass of non-exploding component [Msun]           
- m_2i:  pre-explosion  mass of exploding component     [Msun]       
- a_i:   pre-explosion orbital separation               [cm]
- m_1f:  post-explosion mass of non-exploding component [Msun]           
- m_2f:  post-explosion mass of exploding component     [Msun]   
- vkick: kick velocity                                  [km/s] 
- θ:     polar kick angle (away from e_par)             [rad]
- ϕ:     azimuthal kick angle (off of e_perp)           [rad]
- vimp:  imparted kick velocity on companion            [km/s]     
    
# Output:
- vsys_f: post-explosion systemic velocity            [km/s]
"""
function post_supernova_circular_orbit_vsys(;m_1i, m_2i, a_i, m_1f=-1, m_2f, vkick=0, θ=0, ϕ=0, vimp=0)
    if m_1f == -1
        m_1f = m_1i
    end
    mtilde = (m_1f+m_2f)/(m_1i+m_2i) 
    v_rel = average_orbital_velocity(;m_1=m_1i, m_2=m_2i, a=a_i)
    vkick_div_vrel = vkick/v_rel
    vimp_div_vrel = vimp/v_rel
    # convert trig functions to vars
    cosθ = cos(θ)
    sinθ = sin(θ)
    cosϕ = cos(ϕ)
    sinϕ = sin(ϕ)
    ξ = (1 + vkick_div_vrel^2 + vimp_div_vrel^2 + 2*vkick_div_vrel*cosθ - 2*vimp_div_vrel*sinθ*cosϕ)/(mtilde)
    if (ξ>2)
        return NaN 
    end

    # Systemic velocity
    Δp_x = (m_2f*m_1i - m_2i*m_1f)/(m_2i + m_1i)*vrel + m_2f*vkick*cosθ
    Δp_y = m_1f*vimp + m_2f*vkick*sinθ*cosϕ
    Δp_z = m_2f*vkick*sinθ*sinϕ
    vsys_f = sqrt(Δp_x^2+Δp_y^2+Δp_z^2)/(m_2f+m_1f)/(km)

    return vsys_f
end



###########################################################
### 
### Post kick orbital properties, for general pre-explosion orbits
### 
###########################################################

"""
    post_supernova_general_orbit_parameters(;m_1i, m_2i, a_i, e_i=0, m_1f=-1, m_2f, vkick=0, θ=0, ϕ=0, vimp=0,
        ν_i=0, Ω_i=0, ω_i=0, i_i=0)

Compute post-kick properties for a general pre-explosion orbit 
using equations from [Marchant, Willcox, Vigna-Gomez] TODO

# Arguments:
- m_1i:  pre-explosion  mass of non-exploding component  [Msun]           
- m_2i:  pre-explosion  mass of exploding component      [Msun]       
- a_i:   pre-explosion orbital separation                [cm]
- e_i:   pre-explosion orbital eccentricity              [dimensionless]
- m_1f:  post-explosion mass of non-exploding component  [Msun]           
- m_2f:  post-explosion mass of exploding component      [Msun]   
-
- vkick: kick velocity                                   [km/s] 
- θ:     polar kick angle (away from e_par)              [rad]
- ϕ:     azimuthal kick angle (off of e_perp)            [rad]
- vimp:  imparted kick velocity on companion             [km/s]     
- Orbital orientation angles: 
    - ν_i: true anomaly                                  [rad]
    - Ω_i: pre-explosion longitude of the ascending node [rad]
    - ω_i: pre-explosion argument of periastron          [rad]
    - i_i: pre-explosion inclination                     [rad]
    
# Output: RTW: check! 
- a_f:   post-explosion orbital separation          [cm]
- e_f:   post-explosion orbital eccentricity        [dimensionless]
- Ω_f:   post-explosion longitude of ascending node [rad]      
- ω_f:   post-explosion argument of periastron      [rad]    
- i_f:   post-explosion inclination                 [rad]     
- v_n:   post-explosion systemic velocity, toward N [rad]
- v_e:   post-explosion systemic velocity, toward E [rad]      
- v_rad: post-explosion radial velocity, toward O   [rad]      
"""
function post_supernova_general_orbit_parameters(;m_1i, m_2i, a_i, e_i=0, m_1f=-1, m_2f, vkick=0, 
        θ=0, ϕ=0, vimp=0, ν_i=0, Ω_i=0, ω_i=0, i_i=0)
    if m_1f == -1
        m_1f = m_1i
    end
    M_i = (m_1i+m_2i)*m_sun  
    M_f = (m_1f+m_2f)*m_sun 

    # convert trig functions to vars
    cosθ = cos(θ)
    sinθ = sin(θ)
    cosϕ = cos(ϕ)
    sinϕ = sin(ϕ)
    cosν_i = cos(ν_i)
    sinν_i = sin(ν_i)
    cosi_i = cos(i_i)
    sini_i = sin(i_i)
    cosω_i = cos(ω_i)
    sinω_i = sin(ω_i)
    cosΩ_i = cos(Ω_i)
    sinΩ_i = sin(Ω_i)

    # construct useful intermediary parameters
    f_ν = (1 - e_i^2)/(1 + e_i*cosν_i)
    g_ν = sqrt((1 + 2*e_i*cosν_i + e_i^2)/(1 - e_i^2))

    v_rel = g_ν*sqrt(cgrav*M_i/a_i)/(km)

    h_ν = -e_i*sinν_i/sqrt(1 + 2*e_i*cosν_i + e_i^2)
    j_ν = (1 + e_i*cosν_i)/sqrt(1 + 2*e_i*cosν_i + e_i^2)

    α = vkick/v_rel
    β = vimp/v_rel

    #println("=> ", f_ν, " ", g_ν, " ", h_ν, " ", j_ν, " ")
    #println("==> ", cosθ, " ", cosϕ)
    ξ = f_ν*g_ν^2*M_i/M_f* (1 + α^2 + β^2 + 2* (α*cosθ - h_ν*β* (1 + α*cosθ) - j_ν*β*α*sinθ*cosϕ))
    #println(ξ)

    if ξ>2
        return (NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
    end

    a_f = f_ν*a_i/(2-ξ)
    #println(a_f)

    Lvec_norm = sqrt(α^2*sinθ^2*sinϕ^2+(h_ν*α*sinθ*cosϕ-j_ν*(1+α*cosθ))^2)
    η = f_ν*g_ν^2*M_i/M_f*Lvec_norm^2

    e_f = sqrt(1+(ξ-2)*η)

    #angle between x and direction of motion
    # perhaps can be skipped to save some time, just get directly cos and sin
    τ = acos(max(-1, min(1, j_ν)))
    if (sinν_i < 0)
        τ = -τ
    end
    #δ = ω + ν - τ
    cosτ = cos(τ)
    sinτ = sin(τ)
    cosδ = cosν_i*sinτ*sinω_i - sinν_i*cosτ*sinω_i + sinν_i*sinτ*cosω_i + cosν_i*cosτ*cosω_i
    sinδ = sinν_i*sinτ*sinω_i + cosν_i*cosτ*sinω_i - cosν_i*sinτ*cosω_i + sinν_i*cosτ*cosω_i

    # Elements of rotation matrix
    # RTW: this coordinate system uses E, N, and O, but the text refers instead to celestial W, not E... need to double check
    R_e_par = cosi_i*cosΩ_i*cosδ-sinΩ_i*sinδ
    R_e_per = -(cosi_i*cosΩ_i*sinδ)-sinΩ_i*cosδ
    R_e_z   = -(sini_i*cosΩ_i)
    R_n_par = cosΩ_i*sinδ+cosi_i*sinΩ_i*cosδ
    R_n_per = cosΩ_i*cosδ-cosi_i*sinΩ_i*sinδ
    R_n_z   = -(sini_i*sinΩ_i)
    R_o_par = sini_i*cosδ
    R_o_per = -(sini_i*sinδ)
    R_o_z   = cosi_i

    # velocity, simply compute from change in momentum
    v_par = (-(m_2i-m_2f)*m_1i/M_i*v_rel +(m_1i-m_1f)*m_2i/M_i*v_rel
                    + m_2f*vkick*cosθ + h_ν*m_1f*vimp)/M_f
    v_per = (m_2f*vkick*sinθ*cosϕ+j_ν*m_1f*vimp)/M_f
    v_z = (m_2f*vkick*sinθ*sinϕ)/M_f

    v_e = R_e_par*v_par + R_e_per*v_per + R_e_z*v_z
    v_n = R_n_par*v_par + R_n_per*v_per + R_n_z*v_z
    v_o = R_o_par*v_par + R_o_per*v_per + R_o_z*v_z
    v_rad = -v_o

    # obtain inclination from direction of orbital angular momentum vector
    L_par = j_ν*α*sinθ*sinϕ/Lvec_norm
    L_per = -h_ν*α*sinθ*sinϕ/Lvec_norm
    L_z = h_ν*α*sinθ*cosϕ-j_ν*(1+α*cosθ)/Lvec_norm

    L_e = R_e_par*L_par + R_e_per*L_per + R_e_z*L_z
    L_n = R_n_par*L_par + R_n_per*L_per + R_n_z*L_z
    L_o = R_o_par*L_par + R_o_per*L_per + R_o_z*L_z

    i_f = acos(min(1, abs(L_o)))

    # compute longitude of the ascending node
    # this is obtained by taking the vector n = k x L, where k=(0,0,1) points to the observer
    # then n = (n_e, n_n, 0) points to the ascending node
    n_e = -L_n
    n_n = L_e
    n_norm = sqrt(n_e^2+n_n^2)
    Ω_f = acos(max(-1, min(1, n_n/n_norm)))
    if n_e>0
        Ω_f = 2π - Ω_f
    end
    # compute post-explosion argument of periastron
    if (e_f > 0)
        periastron_angle = acos(max(-1, min(1, 1/e_f*(a_f/(a_i*f_ν)*(1-e_f^2)-1))))
    else
        periastron_angle = 0 # need to check this, though maybe irrelevant as chance of this is null
    end
    # the periastron angle is the same as the true anomaly if the star is moving
    # away from periastron. We need to compute the component of velocity
    # along the line joining both objects (in the COM frame).
    v1y_cm = vimp + h_ν*(m_2i/M_i*v_rel - v_par) - j_ν*v_per
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
        ω_f = angle_to_asc_node-periastron_angle
    else
        ω_f = 2π - (angle_to_asc_node+periastron_angle)
    end
    if ω_f < 0
        ω_f = 2π + ω_f
    end

    return (a_f, e_f, Ω_f, ω_f, i_f, v_n, v_e, v_rad)
end




