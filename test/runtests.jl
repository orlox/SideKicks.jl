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

include("SymbolicFuncs.jl")
create_symbolic_functions_list()

# run a series of tests against the symbolic functions
# the symbolic functions directly compute cross and dot products,
# as well as matrix rotations, so they serve to verify the algebra
# done to obtain the more compact results in the paper

for vkick in [0, 30*km_per_s]
    for vimp in [0, 20*km_per_s]
        for Ω_i in [0.1π, 0.9π, 1.5π, 1.9π]
            for ω_i in [0.1π, 0.9π, 1.5π, 1.9π]
                for i_i in [0.1π/2, 0.5π/2, 0.9π/2]
                    for ν_i in [0.1π, 0.9π, 1.5π, 1.9π]
                        for m1_f in [0.9*m1_i, 0.95*m1_i, m1_i]
                            for m2_f in [0.9*m2_i, 0.95*m2_i, m2_i]
                                a_f_s, e_f_s, v_N_s, v_E_s, v_r_s, Ω_f_s, ω_f_s, i_f_s = 
                                    symbolic_post_kick_parameters_a_e( a=a_i, e=e_i, m_1=m1_i, m_2=m2_i, ν=ν_i, vkick=vkick, 
                                    θ=θ, ϕ=ϕ, v_im=vimp, m_1f = m1_f, m_2f = m2_f, Ω=Ω_i, ω=ω_i, i=i_i, function_list=symbolic_functions_list)
                                a_f, e_f, Ω_f, ω_f, i_f, v_ra, v_dec, v_rad = SideKicks.post_supernova_general_orbit_parameters(m1_i=m1_i, m2_i=m2_i, a_i=a_i, e_i=e_i, m1_f=m1_f, m2_f=m2_f, vkick=vkick, θ=θ, ϕ=ϕ, vimp=vimp, ν_i=ν_i, Ω_i=Ω_i, ω_i=ω_i, i_i=i_i)
                                @test a_f_s ≈ a_f
                                @test e_f_s ≈ e_f
                                @test Ω_f_s ≈ Ω_f atol=1e-6
                                @test ω_f_s ≈ ω_f atol=1e-6
                                @test i_f_s ≈ i_f atol=1e-6
                                @test v_E_s ≈ v_ra atol=1e-6
                                @test v_N_s ≈ v_dec atol=1e-6
                                @test v_r_s ≈ v_rad atol=1e-6

                                if (vimp == 0 && vkick == 0 && m1_f == m1_i && m2_f == m2_i)
                                    # no kick example, verify we recover the input values
                                    @test a_f ≈ a_i
                                    @test e_f ≈ e_i
                                    @test Ω_f ≈ Ω_i atol=1e-6
                                    @test ω_f ≈ ω_i atol=1e-6
                                    @test i_f ≈ i_i atol=1e-6
                                    @test v_ra ≈ 0
                                    @test v_dec ≈ 0
                                    @test v_rad ≈ 0
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
