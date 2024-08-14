using LinearAlgebra
using Symbolics
using IfElse

"""
    symbolic_kick_functions_energy_L()

# Output:
- RTW: TODO!!
"""
function symbolic_kick_functions_energy_L()

    @variables a, e, ν, m1, m1f, m2, m2f, v_par, v_per, v_z, v_im

    M = m1 + m2
    Mf = m1f + m2f

    # semi-major axis of each orbit, star 2 will be the exploding one
    a1 = a*m2/(m1+m2)
    a2 = a*m1/(m1+m2)

    # current vector position with respect to the center of mass for each
    f_ν = (1-e^2)/(1+e*cos(ν))
    r1 = a1*f_ν*[0,1,0]
    r2 = a2*f_ν*[0,-1,0]

    # relative velocity
    g_ν = sqrt((1+2*e*cos(ν)+e^2)/(1-e^2))
    vrel = g_ν*sqrt(cgrav*M/a)
    # vector parallel to motion of exploding star
    e_par = [(1+e*cos(ν)),-e*sin(ν),0]/sqrt(1+2*e*cos(ν)+e^2)
    v1 = -vrel*m2/M*e_par
    v2 = vrel*m1/M*e_par

    #apply kick to star #2
    e_per = [e*sin(ν),(1+e*cos(ν)),0]/sqrt(1+2*e*cos(ν)+e^2) # vector in the orbital plane perpendicular to e_par and pointing into the orbit
    vkick = v_par*e_par + v_per*e_per + v_z*[0,0,1]
    v2f = v2+vkick

    #apply impact to star #1
    v1f = v1 + v_im*[0,1,0]

    #get velocities wrt to center of mass
    v_cm = (m2f*v2f+m1f*v1f)/Mf
    v1cm = v1f-v_cm
    v2cm = v2f-v_cm

    #get position vectors wrt center of mass
    r1cm = a*m2f/Mf*f_ν*[0,1,0]
    r2cm = a*m1f/Mf*f_ν*[0,-1,0]

    #obtain energy
    energy = -cgrav*m1f*m2f/(a*f_ν)+m1f*v1cm⋅v1cm/2+m2f*v2cm⋅v2cm/2
    #obtain angular momentum
    Lvec = m1f*r1cm × v1cm+m2f*r2cm × v2cm

    global energy_function = build_function(energy, [a, e, ν, m1, m1f, m2, m2f, v_par, v_per, v_z, v_im],  expression=Val{false});
    global L_x_function = build_function(Lvec[1],   [a, e, ν, m1, m1f, m2, m2f, v_par, v_per, v_z, v_im],  expression=Val{false});
    global L_y_function = build_function(Lvec[2],   [a, e, ν, m1, m1f, m2, m2f, v_par, v_per, v_z, v_im],  expression=Val{false});
    global L_z_function = build_function(Lvec[3],   [a, e, ν, m1, m1f, m2, m2f, v_par, v_per, v_z, v_im],  expression=Val{false});

    global v_xcm_function = build_function(v_cm[1], [a, e, ν, m1, m1f, m2, m2f, v_par, v_per, v_z, v_im],  expression=Val{false});
    global v_ycm_function = build_function(v_cm[2], [a, e, ν, m1, m1f, m2, m2f, v_par, v_per, v_z, v_im],  expression=Val{false});
    global v_zcm_function = build_function(v_cm[3], [a, e, ν, m1, m1f, m2, m2f, v_par, v_per, v_z, v_im],  expression=Val{false});

    #used to determine if binary moves towards or away from aphelion
    global v_1y_function = build_function(v1cm[2],  [a, e, ν, m1, m1f, m2, m2f, v_par, v_per, v_z, v_im],  expression=Val{false});

    return
end

"""
    symbolic_kick_functions_vcm_and_orbital_elements()

# Arguments:
- RTW: TODO!!

# Output:
- RTW: TODO!!
"""
function symbolic_kick_functions_vcm_and_orbital_elements()
    @variables a, a_final, e, e_final, ν, L_x, L_y, L_z, v_xcm, v_ycm, v_zcm,  v_1y,  Ω,  ω,  i

    f_ν = (1-e^2)/(1+e*cos(ν))

    Rω = [cos(ω+ν) -sin(ω+ν) 0;sin(ω+ν) cos(ω+ν) 0;0 0 1]
    Rι = [cos(i) 0 -sin(i);0 1 0;sin(i) 0 cos(i)]
    RΩ = [cos(Ω) -sin(Ω) 0;sin(Ω) cos(Ω) 0;0 0 1]
    Rtotal  = RΩ*Rι*Rω

    Lvec_rot = Rtotal*[L_x,L_y,L_z]

    Nvec = Lvec_rot×[0,0,1] # vector pointing towards ascending node
    N_norm = sqrt(Nvec⋅Nvec)
    Nvec_norm = Nvec/N_norm 
    Ω_final = IfElse.ifelse(Nvec[1]<=0,acos(Nvec_norm[2]),2*π-acos(Nvec_norm[2]))

    #angle between current position vector and periastron
    #same as true anomaly if binary is moving towards aphelion
    periastron_angle = acos(max(-1,min(1,1/e_final*(a_final/(a*f_ν)*(1-e_final^2)-1))))
    ν_final = IfElse.ifelse(v_1y>=0,periastron_angle,2*π-periastron_angle)

    #apply rotation to unit vector pointing in direction of star #1
    rhat_final = Rtotal*[0,1,0]
    angle_to_asc_node = acos(max(-1,min(1,rhat_final⋅Nvec_norm)))

    y_cross_N_dot_L = (rhat_final × Nvec_norm) ⋅ Lvec_rot

    uncorrected_ω_final = IfElse.ifelse(y_cross_N_dot_L<0,
                                        angle_to_asc_node-ν_final, 2*π-(angle_to_asc_node+ν_final))
    ω_final = IfElse.ifelse(uncorrected_ω_final>0,
        uncorrected_ω_final, 2*π+uncorrected_ω_final)

    #orbital inclination is an easy one
    ι_final = acos(min(1,abs(Lvec_rot[3]/sqrt(L_x^2+L_y^2+L_z^2))))

    #motion of the center of mass
    v_cm_rot = Rtotal*[v_xcm,v_ycm,v_zcm]

    global Ω_function = build_function(Ω_final,  
             [a, a_final, e, e_final, ν, L_x, L_y, L_z, v_xcm, v_ycm, v_zcm, v_1y, Ω, ω, i],  expression=Val{false});
    global ω_function = build_function(ω_final,  
             [a, a_final, e, e_final, ν, L_x, L_y, L_z, v_xcm, v_ycm, v_zcm, v_1y, Ω, ω, i],  expression=Val{false});
    global ι_function = build_function(ι_final, 
             [a, a_final, e, e_final, ν, L_x, L_y, L_z, v_xcm, v_ycm, v_zcm, v_1y, Ω, ω, i],  expression=Val{false});

    global v_N_function = build_function(v_cm_rot[2], 
               [a, a_final, e, e_final, ν, L_x, L_y, L_z, v_xcm, v_ycm, v_zcm, v_1y, Ω, ω, i],  expression=Val{false});
    global v_E_function = build_function(-v_cm_rot[1], 
               [a, a_final, e, e_final, ν, L_x, L_y, L_z, v_xcm, v_ycm, v_zcm, v_1y, Ω, ω, i],  expression=Val{false});
    global v_r_function = build_function(-v_cm_rot[3], 
               [a, a_final, e, e_final, ν, L_x, L_y, L_z, v_xcm, v_ycm, v_zcm, v_1y, Ω, ω, i],  expression=Val{false});

    return
end

"""
    create_symbolic_functions_list()

# Arguments:
- RTW: TODO!!

# Output:
- RTW: TODO!!
"""
function create_symbolic_functions_list()
    symbolic_kick_functions_energy_L()
    symbolic_kick_functions_vcm_and_orbital_elements()
    
    global symbolic_functions_list = (energy_function,L_x_function,L_y_function,L_z_function,v_xcm_function,v_ycm_function,v_zcm_function,v_1y_function,
                            Ω_function, ω_function, ι_function, v_N_function, v_E_function, v_r_function)
    return nothing
end

"""
    wrapped_post_kick_parameters_a_e(a, e, m1, m2, ν, vkick, θ, ϕ, vimp, m1f, m2f, Ω, ω, i)

# Arguments:
- RTW: TODO!!

# Output:
- RTW: TODO!!
"""
function symbolic_post_kick_parameters_a_e(; a, e, m_1, m_2, ν, vkick, θ, ϕ, v_im, m_1f, m_2f, Ω, ω, i, function_list)
    v_par = vkick*cos(θ)
    v_per = vkick*sin(θ)*cos(ϕ)
    v_z = vkick*sin(θ)*sin(ϕ)
    values1 = (a, e, ν, m_1, m_1f, m_2, m_2f, v_par, v_per, v_z, v_im)
    energy = function_list[1](values1) 
    L_x = function_list[2](values1)    
    L_y = function_list[3](values1)    
    L_z = function_list[4](values1)    
    
    v_xcm = function_list[5](values1)  
    v_ycm = function_list[6](values1)  
    v_zcm = function_list[7](values1)  
    v_1y = function_list[8](values1)   

    if energy > 0
        return (NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
    end

    a_final = -cgrav*m_1f*m_2f/(2*energy)
    e_final = sqrt(1-(L_x^2+L_y^2+L_z^2)*(m_1f+m_2f)/(cgrav*a_final*m_1f^2*m_2f^2)+1e-15)

    values2 = (a, a_final, e, e_final, ν, L_x, L_y, L_z, v_xcm, v_ycm, v_zcm, v_1y, Ω, ω, i)

    Ω_f = function_list[9](values2)  
    ω_f = function_list[10](values2) 
    ι_f = function_list[11](values2) 

    v_N = function_list[12](values2) 
    v_E = function_list[13](values2) 
    v_r = function_list[14](values2) 

    return (a_final, e_final, v_N, v_E, v_r, Ω_f, ω_f, ι_f)
end
