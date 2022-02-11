

"""
    kepler_a_from_P(P,m1,m2)

Obtain semimajor axis from period using Kepler's third law

#Arguments:
- P:  orbital period           [d]
- m1: mass of first companion [Msun]
- m2: mass of 2nd companion   [Msun]

#Output:
- a: semi-major axis of the orbit [cm]
"""
function kepler_a_from_P(P,m1,m2)
    return cbrt((P*day_in_sec)^2*cgrav*(m1+m2)*m_sun/(4.0*π^2))
end


"""
    kepler_P_from_a(a,m1,m2)

Obtain period from semimajor axis using Kepler's third law

#Arguments:
- a:  semi-major axis of the orbit [cm]
- m1: mass of first companion     [Msun]
- m2: mass of 2nd companion       [Msun]

#Output:
- P: the orbital period            [d]
"""
function kepler_P_from_a(a,m1,m2)
    return sqrt(4.0*π^2*a^3/(cgrav*(m1+m2)*m_sun))/(day_in_sec)
end


"""
    post_kick_parameters_a(a,m1,m2i,m2f,vkick_kms,theta,phi)

Compute post-kick properties using equations from Tauris et al. (1999)
#Arguments:
- a:     pre-explosion orbital separation            [cm]
- m1:    mass of non-exploding component             [Msun]
- m2i:   mass of exploding component                 [Msun]
- m2f:   post explosion mass of exploding component  [Msun]
- vkick_kms: kick velocity given to m1               [km/s]
- theta: polar angle indicating kick direction       [rad]
- phi:   azimuthal angle indicating kick direction   [rad]
    
#Output:
- afinal: post-explosion orbital separation          [cm]
- efinal: post-explosion excentricity                [dimensionless]
"""
function post_kick_parameters_a(a,mtilde,vkick_div_vrel,theta,phi)
    xi = (1+vkick_div_vrel^2+2*vkick_div_vrel*cos(theta))/(mtilde)
    Q = xi - 1 - (vkick_div_vrel*sin(theta)*cos(phi))^2/(mtilde)

    if (xi>2)
        return (NaN, NaN)
    end
    afinal = a/(2-xi)
    efinal = sqrt(1+(xi-2)*(Q+1))
    return (afinal,efinal)
end


"""
    post_kick_parameters_P(P,m1,m2i,m2f,vkick,theta,phi)

Compute post-kick properties using equations from Tauris et al. (1999), but using the orbital period

#Arguments:
- P:     pre-explosion orbital period               [d]
- m1:    mass of non-exploding component            [Msun]
- m2i:   mass of exploding component                [Msun]
- m2f:   post-explosion mass of exploding component [Msun]
- vkick_kms: kick velocity given to m1              [km/s]
- theta: polar angle indicating kick direction      [rad]
- phi:   azimuthal angle indicating kick direction  [rad]

#Output:
- Ppost: post-explosion orbital period              [d]
- epost: post-explosion eccentricity                [dimensionless]
"""
function post_kick_parameters_P(P,mtilde,vkick_div_vrel,theta,phi)
    xi = (1+vkick_div_vrel^2+2*vkick_div_vrel*cos(theta))/(mtilde)
    Q = xi - 1 - (vkick_div_vrel*sin(theta)*cos(phi))^2/(mtilde)

    if (xi>2)
        return (NaN, NaN)
    end
    Pfinal = P/sqrt(mtilde*(2-xi)^3)
    efinal = sqrt(1+(xi-2)*(Q+1))
    return (Pfinal, efinal)
end

function systemic_velocity(vrel,m1_i,m2_i,m2_f,vkick,θ,ϕ)
    Δp_x = (m2_f*m1_i-m2_i*m1_i)/(m2_i+m1_i)*vrel+m2_f*vkick*cos(θ)
    Δp_y = m2_f*vkick*sin(θ)*cos(ϕ)
    Δp_z = m2_f*vkick*sin(θ)*sin(ϕ)
    return sqrt(Δp_x^2+Δp_y^2+Δp_z^2)/(m2_f+m1_i)/1e5
end

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


