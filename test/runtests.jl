using Test
using SideKicks
using Random
using Turing
include("SymbolicFuncs.jl")

#Test Kepler's laws using earth sun system
year_kepler = SideKicks.kepler_P_from_a(m1=m_sun, m2=m_earth, a=au)
@test year_kepler ≈ year rtol=0.01
au_kepler = SideKicks.kepler_a_from_P(m1=m_sun, m2=m_earth, P=year)
@test au_kepler ≈ au rtol=0.01
earth_vel = SideKicks.relative_velocity(m1=m_sun, m2=m_earth, a=au)
@test earth_vel ≈ 29.7846*km_per_s rtol=0.01

#check K1 and K2 using sun-earth system
RV_s = SideKicks.RV_semiamplitude_K1(m1=m_sun, m2=m_earth, P=year, e=0, i=π/2)
RV_e = SideKicks.RV_semiamplitude_K1(m1=m_earth, m2=m_sun, P=year, e=0, i=π/2)
@test RV_e ≈ 29.7846*km_per_s rtol=0.001 #just test again value that was computed in a trusted version
@test RV_e/RV_s ≈ m_sun/m_earth rtol=0.01 #test that RV ratio follows mass ratio

##########################################################
### 
### Testing circularized orbit kick functions
### 
##########################################################

#Test that a Blaauw kick removing a bit more than half of the mass of a system unbinds the orbit, and one just below does not
orbit_bound =   SideKicks.post_supernova_circular_orbit_a(m1_i=10*m_sun, m2_i=20*m_sun, a_i=au, m2_f=5.0001*m_sun) 
orbit_unbound = SideKicks.post_supernova_circular_orbit_a(m1_i=10*m_sun, m2_i=20*m_sun, a_i=au, m2_f=4.9999*m_sun) 
@test !isnan(orbit_bound[1]) && isnan(orbit_unbound[1])

#Repeat using post_kick_parameters_P
orbit_bound =   SideKicks.post_supernova_circular_orbit_P(m1_i=10*m_sun, m2_i=20*m_sun, P_i=10*day, m2_f=5.0001*m_sun) 
orbit_unbound = SideKicks.post_supernova_circular_orbit_P(m1_i=10*m_sun, m2_i=20*m_sun, P_i=10*day, m2_f=4.9999*m_sun) 
@test !isnan(orbit_bound[1]) && isnan(orbit_unbound[1])

# Generic test values for following kick functions - only use the ones needed
m1_i = 10*m_sun
m2_i = 15*m_sun
m1_f = 0.9*m1_i
m2_f = 0.8*m2_i
vkick=30*km_per_s
vimp=20*km_per_s
θ=π/4
ϕ=π/3
a_i=au
P_i=SideKicks.kepler_P_from_a(m1=m1_i, m2=m2_i, a=a_i)

#Verify consistency between both kick functions when no impact kick
orbit_kick_a = SideKicks.post_supernova_circular_orbit_a(m1_i=m1_i, m2_i=m2_i, a_i=a_i, m2_f=m2_f, vkick=vkick, θ=θ, ϕ=ϕ, vimp=0)
orbit_kick_P = SideKicks.post_supernova_circular_orbit_P(m1_i=m1_i, m2_i=m2_i, P_i=P_i, m2_f=m2_f, vkick=vkick, θ=θ, ϕ=ϕ, vimp=0)
#test period
@test orbit_kick_P[1] ≈ SideKicks.kepler_P_from_a( m1=m1_i, m2=m2_f, a=orbit_kick_a[1]) rtol=0.01
#test eccentricity
@test orbit_kick_a[2]==orbit_kick_P[2]

#Verify consistency between both circular orbit functions (for a and P) with impact kick
orbit_kick_a = SideKicks.post_supernova_circular_orbit_a(m1_i=m1_i, m2_i=m2_i, a_i=a_i, m1_f=m1_f, m2_f=m2_f, vkick=vkick, θ=θ, ϕ=ϕ, vimp=vimp)
orbit_kick_P = SideKicks.post_supernova_circular_orbit_P(m1_i=m1_i, m2_i=m2_i, P_i=P_i, m1_f=m1_f, m2_f=m2_f, vkick=vkick, θ=θ, ϕ=ϕ, vimp=vimp)
#test period
@test orbit_kick_P[1] ≈ SideKicks.kepler_P_from_a( m1=m1_f, m2=m2_f, a=orbit_kick_a[1]) rtol=0.01
#test eccentricity
@test orbit_kick_a[2]==orbit_kick_P[2]



##########################################################
### 
### Testing generalized (non-circular) orbit kick function
### 
##########################################################

#Verify that the generalized kick function produces meaningful output for arbitrary (non disruptive) input
e_i=0.3
ν_i=π/20 
Ω_i=π/20 
ω_i=π/20 
i_i=π/20 
params = SideKicks.post_supernova_general_orbit_parameters(m1_i=m1_i, m2_i=m2_i, a_i=a_i, e_i=e_i, m1_f=m1_f, m2_f=m2_f, vkick=vkick, θ=θ, ϕ=ϕ, vimp=vimp, ν_i=ν_i, Ω_i=Ω_i, ω_i=ω_i, i_i=i_i)
@test all(map(!isnan, params))


#Verify consitency between circularlized and generalized orbit when ecc = 0 and vimp = 0
params_cir_orbit = SideKicks.post_supernova_circular_orbit_a(        m1_i=m1_i, m2_i=m2_i, a_i=a_i,  m2_f=m2_f, vkick=vkick, θ=θ, ϕ=ϕ, vimp=0) 
params_gen_orbit = SideKicks.post_supernova_general_orbit_parameters(m1_i=m1_i, m2_i=m2_i, a_i=a_i,  m2_f=m2_f, vkick=vkick, θ=θ, ϕ=ϕ, vimp=0, e_i=0) 
#test semi-major axis
@test params_cir_orbit[1] ≈ params_gen_orbit[1] rtol=0.01
#test eccentricity
@test params_cir_orbit[2] ≈ params_gen_orbit[2] rtol=0.01

# Check that the generalized kick and the circularized kick are equal when ecc = 0 but vimp != 0
params_cir_orbit = SideKicks.post_supernova_circular_orbit_a(        m1_i=m1_i, m2_i=m2_i, a_i=a_i,  m2_f=m2_f, vkick=vkick, θ=θ, ϕ=ϕ, vimp=vimp) 
params_gen_orbit = SideKicks.post_supernova_general_orbit_parameters(m1_i=m1_i, m2_i=m2_i, a_i=a_i,  m2_f=m2_f, vkick=vkick, θ=θ, ϕ=ϕ, vimp=vimp, e_i=0) 
#test semi-major axis
@test params_cir_orbit[1] ≈ params_gen_orbit[1]
#test eccentricity
@test params_cir_orbit[2] ≈ params_gen_orbit[2]


##########################################################
### 
### Test outputs of all functions against symbolic equivalents
### 
##########################################################

#m1_i = 10
#m2_i = 15
#m2_f = 15
#P_i = 50 # day
#a_i = SideKicks.kepler_a_from_P( m1_i=m1, m2=m2_i, P=P_i)
#vkick = 0
#vim = 0
#θ = 2
#ϕ = 1
#Ω=3
#ω=6
#i=0.4
#e_i=0.4
#ν=2
#
#create_symbolic_functions_list()
#
#a_f, e_f, Ω_f, ω_f, i_f, v_n, v_w, v_rad = 
#    symbolic_post_kick_parameters_a_e( a_i, e_i, m1_i*m_sun, m2_i*m_sun, ν, vkick, 
#    θ, ϕ, vim, m1_i*m_sun, m2_f*m_sun, Ω, ω, i, symbolic_functions_list)
#P_f = SideKicks.kepler_P_from_a( m1=m1, m2=m2_f, a=a_f)
#
#@test isapprox(P, P_f) #test that the period remains the same
#@test isapprox(a, a_f) #test that the semi-major axis remains the same
#@test isapprox(e, e_f, atol=1e-7) #test that the eccentricity remains the same
#@test isapprox(v_N, 0.0, atol=1e-7) #test that the velocity is zero
#@test isapprox(v_E, 0.0, atol=1e-7) #test that the velocity is zero
#@test isapprox(v_r, 0.0, atol=1e-7) #test that the velocity is zero
#@test isapprox(Ω, Ω_f) #test that Ω remains the same
#@test isapprox(ω, ω_f) #test that ω remains the same
#@test isapprox(i, i_f) #test that i remains the same
#
#
##Verify consistency of variables upon no mass loss and no kick
#SideKicks.create_symbolic_functions_list()
#m1 = 10
#m2 = 15
#m2_f = 15
#P = 50 # day
#a = SideKicks.kepler_a_from_P(P,m1,m2)
#vkick = 0
#vim = 0
#θ = 2
#ϕ = 1
#Ω=3
#ω=6
#i=0.4
#e=0.4
#ν=2
#
#a_f, e_f, v_N, v_E, v_r, Ω_f, ω_f, i_f = SideKicks.symbolic_post_kick_parameters_a_e(
#    a,e,m1*m_sun,m2*m_sun,ν,vkick,
#    θ,ϕ,vim, m1*m_sun, m2_f*m_sun,Ω,ω,i,SideKicks.symbolic_functions_list)
#P_f = SideKicks.kepler_P_from_a(a_f,m1,m2_f)
#
#@test isapprox(P, P_f) #test that the period remains the same
#@test isapprox(a, a_f) #test that the semi-major axis remains the same
#@test isapprox(e, e_f, atol=1e-7) #test that the eccentricity remains the same
#@test isapprox(v_N, 0.0, atol=1e-7) #test that the velocity is zero
#@test isapprox(v_E, 0.0, atol=1e-7) #test that the velocity is zero
#@test isapprox(v_r, 0.0, atol=1e-7) #test that the velocity is zero
#@test isapprox(Ω, Ω_f) #test that Ω remains the same
#@test isapprox(ω, ω_f) #test that ω remains the same
#@test isapprox(i, i_f) #test that i remains the same
#
##Verify consistency of variables upon reversing a kick
#SideKicks.create_symbolic_functions_list()
#m1 = 10
#m2 = 15
#m2_f = m2
#P = 50 # day
#a = SideKicks.kepler_a_from_P(P,m1,m2)
#vkick = 100e5
#θ = 2
#ϕ = 1
#Ω=3
#ω=6
#i=0.4
#e=0
#ν=0
#
#a_f1, e_f1, v_N1, v_E1, v_r1, Ω_f1, ω_f1, i_f1 = SideKicks.generalized_post_kick_parameters_a_e(
#    a,e,sin(ν),cos(ν),m1*m_sun,m2*m_sun,m2_f*m_sun,vkick,
#    sin(θ),cos(θ),sin(ϕ),cos(ϕ),sin(Ω),cos(Ω),sin(ω),cos(ω),sin(i),cos(i),SideKicks.symbolic_functions_list)
#
#a_f, e_f, v_N, v_E, v_r, Ω_f, ω_f, i_f = SideKicks.generalized_post_kick_parameters_a_e(
#    a_f1,e_f1,sin(ν),cos(ν),m1*m_sun,m2*m_sun,m2_f*m_sun,-vkick,
#    sin(θ),cos(θ),sin(ϕ),cos(ϕ),sin(Ω_f1),cos(Ω_f1),sin(ω_f1),cos(ω_f1),sin(i_f1),cos(i_f1),SideKicks.symbolic_functions_list)
#P_f = SideKicks.kepler_P_from_a(a_f,m1,m2_f)
#
#@test isapprox(P, P_f, rtol=1e-5) #test that the period remains the same
#@test isapprox(a, a_f, rtol=1e-5) #test that the semi-major axis remains the same
#@test isapprox(e, e_f, atol=1e-7) #test that the eccentricity remains the same
#@test isapprox(v_N1, -v_N, atol=1e-10) #test that the velocity is zero
#@test isapprox(v_E1, -v_E, atol=1e-10) #test that the velocity is zero
#@test isapprox(v_r1, -v_r, atol=1e-10) #test that the velocity is zero
#@test isapprox(Ω, Ω_f) #test that Ω remains the same
#@test isapprox(i, i_f) #test that i remains the same
#
##Verify reversing orbital velocity, with mass loss (Hills 1982)
#SideKicks.create_symbolic_functions_list()
#m1 = 10
#m2 = 15
#m2_f = 10
#P = 50 # day
#a = SideKicks.kepler_a_from_P(P,m1,m2)
#θ = π
#ϕ = 0
#Ω=3
#ω=6
#i=0.4
#e=0.5
#ν=0 #TODO: DOES NOT WORK FOR ALL VALUES OF ν
#r = a*(1-abs2(e))/(1+e*cos(ν))
#vkick = sqrt(cgrav*(m1+m2)*m_sun/a)*2*sqrt(2*a/r-1)
#
#a_f1, e_f1, v_N1, v_E1, v_r1, Ω_f1, ω_f1, i_f1 = SideKicks.symbolic_post_kick_parameters_a_e(
#    a,e,sin(ν),cos(ν),m1*m_sun,m2*m_sun,m2_f*m_sun,0,
#    sin(θ),cos(θ),sin(ϕ),cos(ϕ),sin(Ω),cos(Ω),sin(ω),cos(ω),sin(i),cos(i),SideKicks.symbolic_functions_list)
#P_f1 = SideKicks.kepler_P_from_a(a_f1,m1,m2_f)
#
#a_f2, e_f2, v_N2, v_E2, v_r2, Ω_f2, ω_f2, i_f2 = SideKicks.generalized_post_kick_parameters_a_e(
#    a,e,sin(ν),cos(ν),m1*m_sun,m2*m_sun,m2_f*m_sun,vkick,
#    sin(θ),cos(θ),sin(ϕ),cos(ϕ),sin(Ω),cos(Ω),sin(ω),cos(ω),sin(i),cos(i),SideKicks.symbolic_functions_list)
#P_f2 = SideKicks.kepler_P_from_a(a_f2,m1,m2_f)
#
#@test isapprox(P_f1, P_f2, rtol=1e-5) #test that the period remains the same
#@test isapprox(a_f1, a_f2, rtol=1e-5) #test that the semi-major axis remains the same
#@test isapprox(e_f1, e_f2, atol=1e-7) #test that the eccentricity remains the same
#@test isapprox(Ω_f1+π, Ω_f2) || isapprox(Ω_f1-π, Ω_f2) #test that Ω shifts a factor of π
#@test isapprox(i_f1, i_f2) #test that i remains the same
#
##Verify vanishing orbital velocity, with mass loss (Hills 1982)
#SideKicks.create_symbolic_functions_list()
#m1 = 10
#m2 = 15
#m2_f = 10
#P = 50 # day
#a = SideKicks.kepler_a_from_P(P,m1,m2)
#θ = π
#ϕ = 0
#Ω=3
#ω=6
#i=0.4
#e=0.5
#ν=0 #TODO: DOES NOT WORK FOR ALL VALUES OF ν
#r = a*(1-abs2(e))/(1+e*cos(ν))
#vkick = sqrt(cgrav*(m1+m2)*m_sun/a)*sqrt(2*a/r-1)
#
#a_f, e_f, v_N, v_E, v_r, Ω_f, ω_f, i_f = SideKicks.generalized_post_kick_parameters_a_e(
#    a,e,sin(ν),cos(ν),m1*m_sun,m2*m_sun,m2_f*m_sun,vkick,
#    sin(θ),cos(θ),sin(ϕ),cos(ϕ),sin(Ω),cos(Ω),sin(ω),cos(ω),sin(i),cos(i),SideKicks.symbolic_functions_list)
#P_f = SideKicks.kepler_P_from_a(a_f,m1,m2_f)
#
#@test isapprox(a_f, r/2) #test that the eccentricity becomes 1
#@test isapprox(1., e_f, atol=1e-7) #test that the eccentricity becomes 1
#@test isapprox(Ω, Ω_f) || isapprox(Ω+π, Ω_f) || isapprox(Ω-π, Ω_f) #test that Ω remains the same or flips
#@test isapprox(i, i_f) #test that i remains the same
#
