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


function symbolic_kick_functions_energy_L()

    @variables a,e,m_1,m_2,m_2f,sepν,dνdt,dsepνdν,v_x,v_y,v_z

    # semi-major axis of each orbit, star 2 will be the exploding one
    a_1 = a*m_2/(m_1+m_2)
    a_2 = a*m_1/(m_1+m_2)

    #current vector position with respect to the center of mass for each
    r_1 = a_1*sepν*[0,1,0]#[-sin(ν),cos(ν),0]
    r_2 = a_2*sepν*[0,-1,0]#[sin(ν),-cos(ν),0]

    #velocities
    #v_1 = a_1*dνdt*(dsepνdν*[-sin(ν),cos(ν),0]+sepν*[-cos(ν),-sin(ν),0])
    #v_2 = a_2*dνdt*(dsepνdν*[sin(ν),-cos(ν),0]+sepν*[cos(ν),sin(ν),0])
    v_1 = a_1*dνdt*(dsepνdν*[0,1,0]+sepν*[-1,0,0])
    v_2 = a_2*dνdt*(dsepνdν*[0,-1,0]+sepν*[1,0,0])

    #apply kick to star #2
    #vkick = vk*(sin(θ)*[cos(ν),sin(ν),0]-cos(θ)*cos(ϕ)*[sin(ν),-cos(ν),0]+cos(θ)*cos(ϕ)*[0,0,1])
    vkick = [v_x,v_y,v_z]#vk*(sin(θ)*[1,0,0]-cos(θ)*cos(ϕ)*[0,-1,0]+cos(θ)*cos(ϕ)*[0,0,1])
    v_2f = v_2+vkick

    #get velocities wrt to center of mass
    v_cm = (m_2f*v_2f+m_1*v_1)/(m_2f+m_1)
    v_1cm = v_1-v_cm
    v_2cm = v_2f-v_cm

    #get position vectors wrt center of mass
    r_1cm = a*m_2f/(m_1+m_2f)*sepν*[0,1,0]
    r_2cm = a*m_1/(m_1+m_2f)*sepν*[0,-1,0]

    #obtain energy
    energy = -cgrav*m_1*m_2f/(a*sepν)+m_1*v_1cm⋅v_1cm/2+m_2f*v_2cm⋅v_2cm/2
    #obtain angular momentum
    Lvec = m_1*r_1cm×v_1cm+m_2f*r_2cm×v_2cm

    global energy_function = build_function(energy, [a, e, m_1, m_2, m_2f,sepν,dνdt,dsepνdν, v_x,v_y,v_z], expression=Val{false});
    global L_x_function = build_function(Lvec[1], [a, e, m_1, m_2, m_2f,sepν,dνdt,dsepνdν, v_x,v_y,v_z], expression=Val{false});
    global L_y_function = build_function(Lvec[2], [a, e, m_1, m_2, m_2f,sepν,dνdt,dsepνdν, v_x,v_y,v_z], expression=Val{false});
    global L_z_function = build_function(Lvec[3], [a, e, m_1, m_2, m_2f,sepν,dνdt,dsepνdν, v_x,v_y,v_z], expression=Val{false});

    global v_xcm_function = build_function(v_cm[1], [a, e, m_1, m_2, m_2f,sepν,dνdt,dsepνdν, v_x,v_y,v_z], expression=Val{false});
    global v_ycm_function = build_function(v_cm[2], [a, e, m_1, m_2, m_2f,sepν,dνdt,dsepνdν, v_x,v_y,v_z], expression=Val{false});
    global v_zcm_function = build_function(v_cm[3], [a, e, m_1, m_2, m_2f,sepν,dνdt,dsepνdν, v_x,v_y,v_z], expression=Val{false});

    #used to determine if binary moves towards or away from aphelion
    global v_1y_function = build_function(v_1cm[2], [a, e, m_1, m_2, m_2f,sepν,dνdt,dsepνdν, v_x,v_y,v_z], expression=Val{false});

    return
end

function symbolic_kick_functions_vcm_and_orbital_elements()
    @variables a,a_final,e_final,sepν,L_x,L_y,L_z,v_xcm,v_ycm,v_zcm, v_1y, sinΩ, cosΩ, sinω_plus_ν,cosω_plus_ν, sinι, cosι

    Rω = [cosω_plus_ν -sinω_plus_ν 0;sinω_plus_ν cosω_plus_ν 0;0 0 1]
    Rι = [cosι 0 -sinι;0 1 0;sinι 0 cosι]
    RΩ = [cosΩ -sinΩ 0;sinΩ cosΩ 0;0 0 1]
    Rtotal  = RΩ*Rι*Rω

    Lvec_rot = Rtotal*[L_x,L_y,L_z]

    Nvec = Lvec_rot×[0,0,1]
    N_norm = sqrt(Nvec⋅Nvec)
    Nvec_norm = Nvec/N_norm
    Ω_final = IfElse.ifelse(Nvec[1]<=0,acos(Nvec_norm[2]),2*π-acos(Nvec_norm[2]))

    #angle between current position vector and periastron
    #same as true anomaly if binary is moving towards aphelion
    periastron_angle = acos(max(-1,min(1,1/e_final*(a_final/(a*sepν)*(1-e_final^2)-1))))
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
             [a,a_final,e_final,sepν,L_x,L_y,L_z,v_xcm,v_ycm,v_zcm,v_1y, sinΩ, cosΩ, sinω_plus_ν,cosω_plus_ν, sinι, cosι], expression=Val{false});
    global ω_function = build_function(ω_final, 
             [a,a_final,e_final,sepν,L_x,L_y,L_z,v_xcm,v_ycm,v_zcm,v_1y, sinΩ, cosΩ, sinω_plus_ν,cosω_plus_ν, sinι, cosι], expression=Val{false});
    global ι_function = build_function(ι_final,
             [a,a_final,e_final,sepν,L_x,L_y,L_z,v_xcm,v_ycm,v_zcm,v_1y, sinΩ, cosΩ, sinω_plus_ν,cosω_plus_ν, sinι, cosι], expression=Val{false});

    global v_N_function = build_function(v_cm_rot[2],
               [a,a_final,e_final,sepν,L_x,L_y,L_z,v_xcm,v_ycm,v_zcm,v_1y, sinΩ, cosΩ, sinω_plus_ν,cosω_plus_ν, sinι, cosι], expression=Val{false});
    global v_E_function = build_function(v_cm_rot[1],
               [a,a_final,e_final,sepν,L_x,L_y,L_z,v_xcm,v_ycm,v_zcm,v_1y, sinΩ, cosΩ, sinω_plus_ν,cosω_plus_ν, sinι, cosι], expression=Val{false});
    global v_r_function = build_function(-v_cm_rot[3],
               [a,a_final,e_final,sepν,L_x,L_y,L_z,v_xcm,v_ycm,v_zcm,v_1y, sinΩ, cosΩ, sinω_plus_ν,cosω_plus_ν, sinι, cosι], expression=Val{false});

    return
end

function create_symbolic_functions_list()
    symbolic_kick_functions_energy_L()
    symbolic_kick_functions_vcm_and_orbital_elements()
    
    global symbolic_functions_list = (energy_function,L_x_function,L_y_function,L_z_function,v_xcm_function,v_ycm_function,v_zcm_function,v_1y_function,
                            Ω_function, ω_function, ι_function, v_N_function, v_E_function, v_r_function)
    return nothing
end

function generalized_post_kick_parameters_a_e(a,e,sinν,cosν,m_1,m_2,m_2f,vkick,sinθ,cosθ,sinϕ,cosϕ,sinΩ,cosΩ,sinω,cosω,sinι,cosι,function_list)
    v_x = vkick*cosθ
    v_y = vkick*sinθ*cosϕ
    v_z = vkick*sinθ*sinϕ
    sepν = (1-e^2)/(1+e*cosν)
    dsepνdν = (1-e^2)/(1+e*cosν)^2*(e*sinν)
    dνdt = sqrt(cgrav*(m_1+m_2)/a^3)*(1+e*cosν)^2/sqrt((1-e^2)^3)
    values1 = (a, e, m_1, m_2, m_2f,sepν,dνdt,dsepνdν, v_x,v_y,v_z)
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

    a_final = -cgrav*m_1*m_2f/(2*energy)
    e_final = sqrt(1-(L_x^2+L_y^2+L_z^2)*(m_1+m_2f)/(cgrav*a_final*m_1^2*m_2f^2)+1e-15)

    sinω_plus_ν = sinω*cosν+cosω*sinν
    cosω_plus_ν = cosω*cosν-sinω*sinν
    
    values2 = (a,a_final,e_final,sepν,L_x,L_y,L_z,v_xcm,v_ycm,v_zcm,v_1y, sinΩ,cosΩ, sinω_plus_ν, cosω_plus_ν, sinι, cosι)

    Ω_final = function_list[9](values2) #Ω_function(values2)
    ω_final = function_list[10](values2) #ω_function(values2)
    ι_final = function_list[11](values2) #ι_function(values2)

    v_N = function_list[12](values2) #v_N_function(values2)
    v_E = function_list[13](values2) #v_E_function(values2)
    v_r = function_list[14](values2) #v_r_function(values2)

    return (a_final, e_final, v_N, v_E, v_r, Ω_final, ω_final, ι_final)
end
