using LinearAlgebra
using Symbolics
using IfElse

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
- mtilde: (m1_i+m2_f)/(m1_i+m2_i) [unitless], here we have
    - m1:    mass of non-exploding component
    - m2i:   mass of exploding component
    - m2f:   post explosion mass of exploding component
- vkick_div_vrel: kick velocity in units of the pre-explosion relative velocity    [unitless]
    - vrel = sqrt(G(m1_i + m2_i)/a)
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
- mtilde: (m1_i+m2_f)/(m1_i+m2_i) [unitless], here we have
    - m1:    mass of non-exploding component
    - m2i:   mass of exploding component
    - m2f:   post explosion mass of exploding component
- vkick_div_vrel: kick velocity in units of the pre-explosion relative velocity    [unitless]
    - vrel = sqrt(G(m1_i + m2_i)/a)
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

function post_kick_parameters_a_e(a,e,m1,m2,cosν,sinν,vkick,
    cosθ, sinθ, cosϕ, sinϕ,vim, m1f, m2f,
    cosΩ, sinΩ, cosω, sinω, cosi, sini)

    f_ν = (1-e^2)/(1+e*cosν)
    g_ν = sqrt((1+2*e*cosν+e^2)/(1-e^2))
    M = m1+m2
    Mf = m1f+m2f
    vrel = g_ν*sqrt(cgrav*M/a)

    h_ν = -e*sinν/sqrt(1+2*e*cosν+e^2)
    j_ν = (1+e*cosν)/sqrt(1+2*e*cosν+e^2)

    α = vkick/vrel
    β = vim/vrel

    ξ = f_ν*g_ν^2*M/Mf*(1+α^2+β^2+2*(-h_ν*β*(1+α*cosθ)+α*cosθ-j_ν*β*α*sinθ*cosϕ))

    if ξ>2
        return (NaN, NaN, NaN, NaN, NaN, NaN,NaN,NaN)
    end

    a_f = f_ν*a/(2-ξ)

    Lvec_norm = sqrt(α^2*sinθ^2*sinϕ^2+(h_ν*α*sinθ*cosϕ-j_ν*(1+α*cosθ))^2)
    η = f_ν*g_ν^2*M/Mf*Lvec_norm^2

    e_f = sqrt(1+(ξ-2)*η)

    #angle between x and direction of motion
    # perhaps can be skipped to save some time, just get directly cos and sin
    τ = acos(max(-1,min(1,j_ν)))
    if (sinν < 0)
        τ = -τ
    end
    #δ = ω + ν - τ
    cosτ = cos(τ)
    sinτ = sin(τ)
    cosδ = cosν*sinτ*sinω-sinν*cosτ*sinω+sinν*sinτ*cosω+cosν*cosτ*cosω
    sinδ = sinν*sinτ*sinω+cosν*cosτ*sinω-cosν*sinτ*cosω+sinν*cosτ*cosω

    # Elements of rotation matrix
    R_e_par =cosi*cosΩ*cosδ-sinΩ*sinδ
    R_e_per = -(cosi*cosΩ*sinδ)-sinΩ*cosδ
    R_e_z = -(sini*cosΩ)
    R_n_par = cosΩ*sinδ+cosi*sinΩ*cosδ
    R_n_per = cosΩ*cosδ-cosi*sinΩ*sinδ
    R_n_z = -(sini*sinΩ)
    R_o_par = sini*cosδ
    R_o_per = -(sini*sinδ)
    R_o_z = cosi

    # velocity, simply compute from change in momentum
    v_par = (-(m2-m2f)*m1/M*vrel +(m1-m1f)*m2/M*vrel
                    + m2f*vkick*cosθ + h_ν*m1f*vim)/Mf
    v_per = (m2f*vkick*sinθ*cosϕ+j_ν*m1f*vim)/Mf
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

    i_f = acos(min(1,abs(L_o)))

    # compute longitude of the ascending node
    # this is obtained by taking the vector n = k x L, where k=(0,0,1) points to the observer
    # then n = (n_e, n_n, 0) points to the ascending node
    n_e = -L_n
    n_n = L_e
    n_norm = sqrt(n_e^2+n_n^2)
    Ω_f = acos(max(-1,min(1,n_n/n_norm)))
    if n_e>0
        Ω_f = 2π - Ω_f
    end
    # compute post-explosion argument of periastron
    if (e_f > 0)
        periastron_angle = acos(max(-1,min(1,1/e_f*(a_f/(a*f_ν)*(1-e_f^2)-1))))
    else
        periastron_angle = 0 # need to check this, though maybe irrelevant as chance of this is null
    end
    # the periastron angle is the same as the true anomaly if the star is moving
    # away from periastron. We need to compute the component of velocity
    # amont the line joining both objects (in the COM frame).
    v1y_cm = h_ν*(m2/M*vrel - v_par) - j_ν*v_per
    if v1y_cm>0
        periastron_angle = 2π-periastron_angle
    end
    # compute the angle from current location to ascending node
    # in the rotated frame. This is R*\hat{y}
    rvec_e = R_e_par*h_ν + R_e_per*j_ν
    rvec_n = R_n_par*h_ν + R_n_per*j_ν
    rvec_o = R_o_par*h_ν + R_o_per*j_ν
    dot_prod = (n_e*rvec_e + n_n*rvec_n)/n_norm
    angle_to_asc_node = acos(max(-1,min(1,dot_prod)))
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

    return (a_f, e_f, v_n, v_e, -v_o, Ω_f,ω_f,i_f)
end

function wrapped_post_kick_parameters_a_e(a,e,m1,m2,ν,vkick, θ, ϕ,vim, m1f, m2f, Ω, ω, i)
    cosν = cos(ν)
    sinν = sin(ν)
    cosθ = cos(θ)
    sinθ = sin(θ)
    cosϕ = cos(ϕ)
    sinϕ = sin(ϕ)

    cosi = cos(i)
    sini = sin(i)
    cosω = cos(ω)
    sinω = sin(ω)
    cosΩ = cos(Ω)
    sinΩ = sin(Ω)

    return post_kick_parameters_a_e(a,e,m1,m2,cosν,sinν,vkick,
        cosθ, sinθ, cosϕ, sinϕ,vim, m1f, m2f,
        cosΩ, sinΩ, cosω, sinω, cosi, sini)
end

function symbolic_kick_functions_energy_L()

    @variables a,e,ν,m_1,m_1f,m_2,m_2f,v_par,v_per,v_z,v_im

    M = m_1 + m_2
    Mf = m_1f + m_2f

    # semi-major axis of each orbit, star 2 will be the exploding one
    a_1 = a*m_2/(m_1+m_2)
    a_2 = a*m_1/(m_1+m_2)

    # current vector position with respect to the center of mass for each
    f_ν = (1-e^2)/(1+e*cos(ν))
    r_1 = a_1*f_ν*[0,1,0]
    r_2 = a_2*f_ν*[0,-1,0]

    # relative velocity
    g_ν = sqrt((1+2*e*cos(ν)+e^2)/(1-e^2))
    vrel = g_ν*sqrt(cgrav*M/a)
    # vector parallel to motion of exploding star
    e_par = [(1+e*cos(ν)),-e*sin(ν),0]/sqrt(1+2*e*cos(ν)+e^2)
    v_1 = -vrel*m_2/M*e_par
    v_2 = vrel*m_1/M*e_par

    #apply kick to star #2
    e_per = [e*sin(ν),(1+e*cos(ν)),0]/sqrt(1+2*e*cos(ν)+e^2) # vector in the orbital plane perpendicular to e_par and pointing into the orbit
    vkick = v_par*e_par + v_per*e_per + v_z*[0,0,1]
    v_2f = v_2+vkick

    #apply impact to star #1
    v_1f = v_1 + v_im*[0,1,0]

    #get velocities wrt to center of mass
    v_cm = (m_2f*v_2f+m_1f*v_1f)/Mf
    v_1cm = v_1f-v_cm
    v_2cm = v_2f-v_cm

    #get position vectors wrt center of mass
    r_1cm = a*m_2f/Mf*f_ν*[0,1,0]
    r_2cm = a*m_1f/Mf*f_ν*[0,-1,0]

    #obtain energy
    energy = -cgrav*m_1f*m_2f/(a*f_ν)+m_1f*v_1cm⋅v_1cm/2+m_2f*v_2cm⋅v_2cm/2
    #obtain angular momentum
    Lvec = m_1f*r_1cm×v_1cm+m_2f*r_2cm×v_2cm

    global energy_function = build_function(energy, [a,e,ν,m_1,m_1f,m_2,m_2f,v_par,v_per,v_z,v_im], expression=Val{false});
    global L_x_function = build_function(Lvec[1], [a,e,ν,m_1,m_1f,m_2,m_2f,v_par,v_per,v_z,v_im], expression=Val{false});
    global L_y_function = build_function(Lvec[2], [a,e,ν,m_1,m_1f,m_2,m_2f,v_par,v_per,v_z,v_im], expression=Val{false});
    global L_z_function = build_function(Lvec[3], [a,e,ν,m_1,m_1f,m_2,m_2f,v_par,v_per,v_z,v_im], expression=Val{false});

    global v_xcm_function = build_function(v_cm[1], [a,e,ν,m_1,m_1f,m_2,m_2f,v_par,v_per,v_z,v_im], expression=Val{false});
    global v_ycm_function = build_function(v_cm[2], [a,e,ν,m_1,m_1f,m_2,m_2f,v_par,v_per,v_z,v_im], expression=Val{false});
    global v_zcm_function = build_function(v_cm[3], [a,e,ν,m_1,m_1f,m_2,m_2f,v_par,v_per,v_z,v_im], expression=Val{false});

    #used to determine if binary moves towards or away from aphelion
    global v_1y_function = build_function(v_1cm[2], [a,e,ν,m_1,m_1f,m_2,m_2f,v_par,v_per,v_z,v_im], expression=Val{false});

    return
end

function symbolic_kick_functions_vcm_and_orbital_elements()
    @variables a,a_final,e,e_final,ν,L_x,L_y,L_z,v_xcm,v_ycm,v_zcm, v_1y, Ω, ω, i

    f_ν = (1-e^2)/(1+e*cos(ν))

    Rω = [cos(ω+ν) -sin(ω+ν) 0;sin(ω+ν) cos(ω+ν) 0;0 0 1]
    Rι = [cos(i) 0 -sin(i);0 1 0;sin(i) 0 cos(i)]
    RΩ = [cos(Ω) -sin(Ω) 0;sin(Ω) cos(Ω) 0;0 0 1]
    Rtotal  = RΩ*Rι*Rω

    Lvec_rot = Rtotal*[L_x,L_y,L_z]

    Nvec = Lvec_rot×[0,0,1]
    N_norm = sqrt(Nvec⋅Nvec)
    Nvec_norm = Nvec/N_norm
    Ω_final = IfElse.ifelse(Nvec[1]<=0,acos(Nvec_norm[2]),2*π-acos(Nvec_norm[2]))

    #angle between current position vector and periastron
    #same as true anomaly if binary is moving towards aphelion
    periastron_angle = acos(max(-1,min(1,1/e_final*(a_final/(a*f_ν)*(1-e_final^2)-1))))
    ν_final = IfElse.ifelse(v_1y>=0,periastron_angle,2*π-periastron_angle)

    #apply rotation to unit vector pointing in direction of star #1
    rhat_final = Rtotal*[0,1,0]
    L_cross_L_cross_N = Lvec_rot×(Lvec_rot×[0,0,1])
    angle_to_asc_node = acos(max(-1,min(1,rhat_final⋅Nvec_norm)))
    uncorrected_ω_final = IfElse.ifelse(L_cross_L_cross_N⋅rhat_final>0,
                                        angle_to_asc_node-ν_final, 2*π-(angle_to_asc_node+ν_final))
    ω_final = IfElse.ifelse(uncorrected_ω_final>0,
        uncorrected_ω_final, 2*π+uncorrected_ω_final)

    #orbital inclination is an easy one
    ι_final = acos(min(1,abs(Lvec_rot[3]/sqrt(L_x^2+L_y^2+L_z^2))))

    #motion of the center of mass
    v_cm_rot = Rtotal*[v_xcm,v_ycm,v_zcm]

    global Ω_function = build_function(Ω_final, 
             [a,a_final,e,e_final,ν,L_x,L_y,L_z,v_xcm,v_ycm,v_zcm,v_1y, Ω, ω, i], expression=Val{false});
    global ω_function = build_function(ω_final, 
             [a,a_final,e,e_final,ν,L_x,L_y,L_z,v_xcm,v_ycm,v_zcm,v_1y, Ω, ω, i], expression=Val{false});
    global ι_function = build_function(ι_final,
             [a,a_final,e,e_final,ν,L_x,L_y,L_z,v_xcm,v_ycm,v_zcm,v_1y, Ω, ω, i], expression=Val{false});

    global v_N_function = build_function(v_cm_rot[2],
               [a,a_final,e,e_final,ν,L_x,L_y,L_z,v_xcm,v_ycm,v_zcm,v_1y, Ω, ω, i], expression=Val{false});
    global v_E_function = build_function(v_cm_rot[1],
               [a,a_final,e,e_final,ν,L_x,L_y,L_z,v_xcm,v_ycm,v_zcm,v_1y, Ω, ω, i], expression=Val{false});
    global v_r_function = build_function(-v_cm_rot[3],
               [a,a_final,e,e_final,ν,L_x,L_y,L_z,v_xcm,v_ycm,v_zcm,v_1y, Ω, ω, i], expression=Val{false});

    return
end

function create_symbolic_functions_list()
    symbolic_kick_functions_energy_L()
    symbolic_kick_functions_vcm_and_orbital_elements()
    
    global symbolic_functions_list = (energy_function,L_x_function,L_y_function,L_z_function,v_xcm_function,v_ycm_function,v_zcm_function,v_1y_function,
                            Ω_function, ω_function, ι_function, v_N_function, v_E_function, v_r_function)
    return nothing
end
function symbolic_post_kick_parameters_a_e(a,e,m_1,m_2,ν,vkick, θ, ϕ,v_im, m_1f, m_2f, Ω, ω, i, function_list)
    v_par = vkick*cos(θ)
    v_per = vkick*sin(θ)*cos(ϕ)
    v_z = vkick*sin(θ)*sin(ϕ)
    values1 = (a,e,ν,m_1,m_1f,m_2,m_2f,v_par,v_per,v_z,v_im)
    energy = function_list[1](values1) #energy_function(valuesa)
    L_x = function_list[2](values1) #L_x_function(valuesa)
    L_y = function_list[3](values1) #L_y_function(valuesa)
    L_z = function_list[4](values1) #L_z_function(valuesa)
    
    v_xcm = function_list[5](values1) #v_xcm_function(valuesa)
    v_ycm = function_list[6](values1) #v_ycm_function(valuesa)
    v_zcm = function_list[7](values1) #v_zcm_function(valuesa)
    v_1y = function_list[8](values1) #v_1y_function(valuesa)

    if energy > 0
        return (NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
    end

    a_final = -cgrav*m_1f*m_2f/(2*energy)
    e_final = sqrt(1-(L_x^2+L_y^2+L_z^2)*(m_1f+m_2f)/(cgrav*a_final*m_1f^2*m_2f^2)+1e-15)

    values2 = (a,a_final,e,e_final,ν,L_x,L_y,L_z,v_xcm,v_ycm,v_zcm,v_1y, Ω, ω, i)

    Ω_final = function_list[9](values2) #Ω_function(values2)
    ω_final = function_list[10](values2) #ω_function(values2)
    ι_final = function_list[11](values2) #ι_function(values2)

    v_N = function_list[12](values2) #v_N_function(values2)
    v_E = function_list[13](values2) #v_E_function(values2)
    v_r = function_list[14](values2) #v_r_function(values2)

    return (a_final, e_final, v_N, v_E, v_r, Ω_final, ω_final, ι_final)
end
