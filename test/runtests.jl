using Test
using SideKicks
using Random
using Turing
include("SymbolicFuncs.jl")

#Test Kepler's laws using earth sun system
year_in_days = SideKicks.kepler_P_from_a(;m_1=1, m_2=m_earth/m_sun, a=au)
@test abs(year_in_days-365.25) < 0.01
au_from_kepler = SideKicks.kepler_a_from_P(;m_1=1, m_2=m_earth/m_sun, P=year_in_days)
@test isapprox(au, au_from_kepler)
earth_vel = SideKicks.average_orbital_velocity(;m_1=1, m_2=m_earth/m_sun, a=au)
@test abs(earth_vel-29.7846) < 0.01

#check K1 and K2 using sun-earth system
RV_e = SideKicks.RV_semiamplitude_K1(;m_1=m_earth/m_sun, m_2=1, P=year_in_days, e=0, i=π/2)
RV_s = SideKicks.RV_semiamplitude_K1(;m_1=1, m_2=m_earth/m_sun, P=year_in_days, e=0, i=π/2)
@test abs(RV_e-29.7846) < 0.001 #just test agains value that was computed in a trusted version
@test isapprox(RV_e/RV_s, m_sun/m_earth) #test that RV ratio follows mass ratio

##########################################################
### 
### Testing circularized orbit kick functions
### 
##########################################################

#Test that a Blaauw kick removing a bit more than half of the mass of a system unbinds the orbit, and one just below does not
orbit_bound =   SideKicks.post_supernova_circular_orbit_a(;m_1i=10, m_2i=20, a_i=au, m_2f=5.0001) 
orbit_unbound = SideKicks.post_supernova_circular_orbit_a(;m_1i=10, m_2i=20, a_i=au, m_2f=4.9999) 
@test !isnan(orbit_bound[1]) && isnan(orbit_unbound[1])

#Repeat using post_kick_parameters_P
orbit_bound =   SideKicks.post_supernova_circular_orbit_P(;m_1i=10, m_2i=20, P_i=10, m_2f=5.0001) 
orbit_unbound = SideKicks.post_supernova_circular_orbit_P(;m_1i=10, m_2i=20, P_i=10, m_2f=4.9999) 
@test !isnan(orbit_bound[1]) && isnan(orbit_unbound[1])

# Generic test values for following kick functions - only use the ones needed
m_1i = 10
m_2i = 15
m_1f = 0.9*m_1i
m_2f = 0.8*m_2i
vkick=30 
vimp=20
θ=π/4
ϕ=π/3
a_i=au
P_i=SideKicks.kepler_P_from_a(;m_1=m_1i, m_2=m_2i, a=a_i)

#Verify consistency between both kick functions when no impact kick
orbit_kick_a = SideKicks.post_supernova_circular_orbit_a(;m_1i=m_1i, m_2i=m_2i, a_i=a_i, m_2f=m_2f, vkick=vkick, θ=θ, ϕ=ϕ, vimp=0)
orbit_kick_P = SideKicks.post_supernova_circular_orbit_P(;m_1i=m_1i, m_2i=m_2i, P_i=P_i, m_2f=m_2f, vkick=vkick, θ=θ, ϕ=ϕ, vimp=0)
#test period
@test isapprox(orbit_kick_P[1], SideKicks.kepler_P_from_a(; m_1=m_1i, m_2=m_2f, a=orbit_kick_a[1]))
#test eccentricity
@test orbit_kick_a[2]==orbit_kick_P[2]

#Verify consistency between both circular orbit functions (for a and P) with impact kick
orbit_kick_a = SideKicks.post_supernova_circular_orbit_a(;m_1i=m_1i, m_2i=m_2i, a_i=a_i, m_1f=m_1f, m_2f=m_2f, vkick=vkick, θ=θ, ϕ=ϕ, vimp=vimp)
orbit_kick_P = SideKicks.post_supernova_circular_orbit_P(;m_1i=m_1i, m_2i=m_2i, P_i=P_i, m_1f=m_1f, m_2f=m_2f, vkick=vkick, θ=θ, ϕ=ϕ, vimp=vimp)
#test period
@test isapprox(orbit_kick_P[1], SideKicks.kepler_P_from_a(; m_1=m_1f, m_2=m_2f, a=orbit_kick_a[1]))
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
params = SideKicks.post_supernova_general_orbit_parameters(;m_1i=m_1i, m_2i=m_2i, a_i=a_i, e_i=e_i, m_1f=m_1f, m_2f=m_2f, vkick=vkick, θ=θ, ϕ=ϕ, vimp=vimp, ν_i=ν_i, Ω_i=Ω_i, ω_i=ω_i, i_i=i_i)
@test all(map(!isnan, params))


#Verify consitency between circularlized and generalized orbit when ecc = 0 and vimp = 0
params_cir_orbit = SideKicks.post_supernova_circular_orbit_a(        ;m_1i=m_1i, m_2i=m_2i, a_i=a_i,  m_2f=m_2f, vkick=vkick, θ=θ, ϕ=ϕ, vimp=0) 
params_gen_orbit = SideKicks.post_supernova_general_orbit_parameters(;m_1i=m_1i, m_2i=m_2i, a_i=a_i,  m_2f=m_2f, vkick=vkick, θ=θ, ϕ=ϕ, vimp=0, e_i=0) 
#test semi-major axis
@test isapprox(params_cir_orbit[1], params_gen_orbit[1])
#test eccentricity
@test isapprox(params_cir_orbit[2], params_gen_orbit[2])

# Check that the generalized kick and the circularized kick are equal when ecc = 0 but vimp != 0
params_cir_orbit = SideKicks.post_supernova_circular_orbit_a(        ;m_1i=m_1i, m_2i=m_2i, a_i=a_i,  m_2f=m_2f, vkick=vkick, θ=θ, ϕ=ϕ, vimp=vimp) 
params_gen_orbit = SideKicks.post_supernova_general_orbit_parameters(;m_1i=m_1i, m_2i=m_2i, a_i=a_i,  m_2f=m_2f, vkick=vkick, θ=θ, ϕ=ϕ, vimp=vimp, e_i=0) 
#test semi-major axis
@test isapprox(params_cir_orbit[1], params_gen_orbit[1])
#test eccentricity
@test isapprox(params_cir_orbit[2], params_gen_orbit[2])


##########################################################
### 
### Test outputs of all functions against symbolic equivalents
### 
##########################################################

m1 = 10
m2 = 15
m2f = 15
P_i = 50 # days
a_i = SideKicks.kepler_a_from_P( ;m_1=m1, m_2=m2, P=P_i)
vkick = 0
vim = 0
θ = 2
ϕ = 1
Ω_i=3
ω_i=6
i_i=0.4
e_i=0.4
ν=2

create_symbolic_functions_list()

a_f, e_f, v_N, v_E, v_r, Ω_f, ω_f, i_f = 
    symbolic_post_kick_parameters_a_e( a_i, e_i, m1*m_sun, m2*m_sun, ν, vkick, 
    θ, ϕ, vim, m1*m_sun, m2f*m_sun, Ω_i, ω_i, i_i, symbolic_functions_list)
P_f = SideKicks.kepler_P_from_a( ;m_1=m1, m_2=m2f, a=a_f)

@test isapprox(P_i, P_f) #test that the period remains the same
@test isapprox(a_i, a_f) #test that the semi-major axis remains the same
@test isapprox(e_i, e_f, atol=1e-7) #test that the eccentricity remains the same
@test isapprox(v_N, 0.0, atol=1e-7) #test that the velocity is zero
@test isapprox(v_E, 0.0, atol=1e-7) #test that the velocity is zero
@test isapprox(v_r, 0.0, atol=1e-7) #test that the velocity is zero
@test isapprox(Ω_i, Ω_f) #test that Ω remains the same
@test isapprox(ω_i, ω_f) #test that ω remains the same
@test isapprox(i_i, i_f) #test that ι remains the same


#Verify consistency of variables upon no mass loss and no kick
SideKicks.create_symbolic_functions_list()
m1 = 10
m2 = 15
m2f = 15
P_i = 50 # days
a_i = SideKicks.kepler_a_from_P(P_i,m1,m2)
vkick = 0
vim = 0
θ = 2
ϕ = 1
Ω_i=3
ω_i=6
i_i=0.4
e_i=0.4
ν=2

a_f, e_f, v_N, v_E, v_r, Ω_f, ω_f, i_f = SideKicks.symbolic_post_kick_parameters_a_e(
    a_i,e_i,m1*m_sun,m2*m_sun,ν,vkick,
    θ,ϕ,vim, m1*m_sun, m2f*m_sun,Ω_i,ω_i,i_i,SideKicks.symbolic_functions_list)
P_f = SideKicks.kepler_P_from_a(a_f,m1,m2f)

@test isapprox(P_i, P_f) #test that the period remains the same
@test isapprox(a_i, a_f) #test that the semi-major axis remains the same
@test isapprox(e_i, e_f, atol=1e-7) #test that the eccentricity remains the same
@test isapprox(v_N, 0.0, atol=1e-7) #test that the velocity is zero
@test isapprox(v_E, 0.0, atol=1e-7) #test that the velocity is zero
@test isapprox(v_r, 0.0, atol=1e-7) #test that the velocity is zero
@test isapprox(Ω_i, Ω_f) #test that Ω remains the same
@test isapprox(ω_i, ω_f) #test that ω remains the same
@test isapprox(i_i, i_f) #test that ι remains the same

#Verify consistency of variables upon reversing a kick
SideKicks.create_symbolic_functions_list()
m1 = 10
m2 = 15
m2f = m2
P_i = 50 # days
a_i = SideKicks.kepler_a_from_P(P_i,m1,m2)
vkick = 100e5
θ = 2
ϕ = 1
Ω_i=3
ω_i=6
ι_i=0.4
e_i=0
ν=0

a_f1, e_f1, v_N1, v_E1, v_r1, Ω_f1, ω_f1, ι_f1 = SideKicks.generalized_post_kick_parameters_a_e(
    a_i,e_i,sin(ν),cos(ν),m1*m_sun,m2*m_sun,m2f*m_sun,vkick,
    sin(θ),cos(θ),sin(ϕ),cos(ϕ),sin(Ω_i),cos(Ω_i),sin(ω_i),cos(ω_i),sin(ι_i),cos(ι_i),SideKicks.symbolic_functions_list)

a_f, e_f, v_N, v_E, v_r, Ω_f, ω_f, ι_f = SideKicks.generalized_post_kick_parameters_a_e(
    a_f1,e_f1,sin(ν),cos(ν),m1*m_sun,m2*m_sun,m2f*m_sun,-vkick,
    sin(θ),cos(θ),sin(ϕ),cos(ϕ),sin(Ω_f1),cos(Ω_f1),sin(ω_f1),cos(ω_f1),sin(ι_f1),cos(ι_f1),SideKicks.symbolic_functions_list)
P_f = SideKicks.kepler_P_from_a(a_f,m1,m2f)

@test isapprox(P_i, P_f, rtol=1e-5) #test that the period remains the same
@test isapprox(a_i, a_f, rtol=1e-5) #test that the semi-major axis remains the same
@test isapprox(e_i, e_f, atol=1e-7) #test that the eccentricity remains the same
@test isapprox(v_N1, -v_N, atol=1e-10) #test that the velocity is zero
@test isapprox(v_E1, -v_E, atol=1e-10) #test that the velocity is zero
@test isapprox(v_r1, -v_r, atol=1e-10) #test that the velocity is zero
@test isapprox(Ω_i, Ω_f) #test that Ω remains the same
@test isapprox(ι_i, ι_f) #test that ι remains the same

#Verify reversing orbital velocity, with mass loss (Hills 1982)
SideKicks.create_symbolic_functions_list()
m1 = 10
m2 = 15
m2f = 10
P_i = 50 # days
a_i = SideKicks.kepler_a_from_P(P_i,m1,m2)
θ = π
ϕ = 0
Ω_i=3
ω_i=6
ι_i=0.4
e_i=0.5
ν=0 #TODO: DOES NOT WORK FOR ALL VALUES OF ν
r_i = a_i*(1-abs2(e_i))/(1+e_i*cos(ν))
vkick = sqrt(cgrav*(m1+m2)*m_sun/a_i)*2*sqrt(2*a_i/r_i-1)

a_f1, e_f1, v_N1, v_E1, v_r1, Ω_f1, ω_f1, ι_f1 = SideKicks.symbolic_post_kick_parameters_a_e(
    a_i,e_i,sin(ν),cos(ν),m1*m_sun,m2*m_sun,m2f*m_sun,0,
    sin(θ),cos(θ),sin(ϕ),cos(ϕ),sin(Ω_i),cos(Ω_i),sin(ω_i),cos(ω_i),sin(ι_i),cos(ι_i),SideKicks.symbolic_functions_list)
P_f1 = SideKicks.kepler_P_from_a(a_f1,m1,m2f)

a_f2, e_f2, v_N2, v_E2, v_r2, Ω_f2, ω_f2, ι_f2 = SideKicks.generalized_post_kick_parameters_a_e(
    a_i,e_i,sin(ν),cos(ν),m1*m_sun,m2*m_sun,m2f*m_sun,vkick,
    sin(θ),cos(θ),sin(ϕ),cos(ϕ),sin(Ω_i),cos(Ω_i),sin(ω_i),cos(ω_i),sin(ι_i),cos(ι_i),SideKicks.symbolic_functions_list)
P_f2 = SideKicks.kepler_P_from_a(a_f2,m1,m2f)

@test isapprox(P_f1, P_f2, rtol=1e-5) #test that the period remains the same
@test isapprox(a_f1, a_f2, rtol=1e-5) #test that the semi-major axis remains the same
@test isapprox(e_f1, e_f2, atol=1e-7) #test that the eccentricity remains the same
@test isapprox(Ω_f1+π, Ω_f2) || isapprox(Ω_f1-π, Ω_f2) #test that Ω shifts a factor of π
@test isapprox(ι_f1, ι_f2) #test that ι remains the same

#Verify vanishing orbital velocity, with mass loss (Hills 1982)
SideKicks.create_symbolic_functions_list()
m1 = 10
m2 = 15
m2f = 10
P_i = 50 # days
a_i = SideKicks.kepler_a_from_P(P_i,m1,m2)
θ = π
ϕ = 0
Ω_i=3
ω_i=6
ι_i=0.4
e_i=0.5
ν=0 #TODO: DOES NOT WORK FOR ALL VALUES OF ν
r_i = a_i*(1-abs2(e_i))/(1+e_i*cos(ν))
vkick = sqrt(cgrav*(m1+m2)*m_sun/a_i)*sqrt(2*a_i/r_i-1)

a_f, e_f, v_N, v_E, v_r, Ω_f, ω_f, ι_f = SideKicks.generalized_post_kick_parameters_a_e(
    a_i,e_i,sin(ν),cos(ν),m1*m_sun,m2*m_sun,m2f*m_sun,vkick,
    sin(θ),cos(θ),sin(ϕ),cos(ϕ),sin(Ω_i),cos(Ω_i),sin(ω_i),cos(ω_i),sin(ι_i),cos(ι_i),SideKicks.symbolic_functions_list)
P_f = SideKicks.kepler_P_from_a(a_f,m1,m2f)

@test isapprox(a_f, r_i/2) #test that the eccentricity becomes 1
@test isapprox(1., e_f, atol=1e-7) #test that the eccentricity becomes 1
@test isapprox(Ω_i, Ω_f) || isapprox(Ω_i+π, Ω_f) || isapprox(Ω_i-π, Ω_f) #test that Ω remains the same or flips
@test isapprox(ι_i, ι_f) #test that ι remains the same

