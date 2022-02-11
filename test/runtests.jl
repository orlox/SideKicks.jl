using Test
using SideKicks
using Random
using Turing

#Test Kepler's laws using earth sun system
year_in_days = SideKicks.kepler_P_from_a(au,1,m_earth/m_sun)
@test abs(year_in_days-365.25) < 0.01
au_from_kepler = SideKicks.kepler_a_from_P(year_in_days,1,m_earth/m_sun)
@test isapprox(au, au_from_kepler)

#Test that a Blaauw kick removing a bit more than half of the mass of a system unbinds the orbit, and one just below does not
orbit_bound = SideKicks.post_kick_parameters_P(10,0.4999,0,0,0)
orbit_unbound = SideKicks.post_kick_parameters_P(10,0.5001,0,0,0)
@test isnan(orbit_bound[1]) && !isnan(orbit_unbound[1])

#Repeat using post_kick_parameters_a
orbit_bound = SideKicks.post_kick_parameters_a(au,0.4999,0,0,0)
orbit_unbound = SideKicks.post_kick_parameters_a(au,0.5001,0,0,0)
@test isnan(orbit_bound[1]) && !isnan(orbit_unbound[1])

#Verify consistency between both kick functions
m1_i = 10
m2_i = 15
frac = 0.8
m2_f = m2_i*frac
mtilde = (m1_i+m2_f)/(m1_i+m2_i)
vkick_div_vrel = 0.1
θ=π/4
ϕ=π/3
semimajor_axis = au
period_in_days = SideKicks.kepler_P_from_a(au,m1_i,m2_i)
orbit_kick_a = SideKicks.post_kick_parameters_a(au,mtilde,vkick_div_vrel,θ,ϕ)
orbit_kick_P = SideKicks.post_kick_parameters_P(period_in_days,mtilde,vkick_div_vrel,θ,ϕ)
#test period
@test isapprox(SideKicks.kepler_P_from_a(orbit_kick_a[1],m1_i,m2_f), orbit_kick_P[1])
#test eccentricity
@test orbit_kick_a[2]==orbit_kick_P[2]

#check K1 and K2 using sun-earth system
RV_e = SideKicks.RV_semiamplitude_K(year_in_days,0,π/2,m_earth/m_sun,1)
RV_s = SideKicks.RV_semiamplitude_K(year_in_days,0,π/2,1,m_earth/m_sun)
@test abs(RV_e-29.7846) < 0.001 #just test agains value that was computed in a trusted version
@test isapprox(RV_e/RV_s, m_sun/m_earth) #test that RV ratio follows mass ratio

#do sampling on a simplified problem to see if results are consistent
#use binary system tested above with narrow priors on kick direction and fraction of
#mass loss to make problem invertible
Random.seed!(1234);
period = orbit_kick_P[1]
eccentricity = orbit_kick_P[2]
K1 = SideKicks.RV_semiamplitude_K(period,eccentricity,π/3,m1_i,m2_f)
K2 = SideKicks.RV_semiamplitude_K(period,eccentricity,π/3,m2_f,m1_i)
test_model = SideKicks.createMCMCModel([:P,:e,:K1,:K2],
                         [period,eccentricity,K1,K2],
                         [period,eccentricity,K1,K2]*0.001,
                        logm1_i_dist = Normal(log10(m1_i),0.001),
                        frac_dist = truncated(Normal(frac,0.001),0.0001,0.9999),
                        cos_θ_dist = truncated(Normal(cos(θ),0.001),-1,1),
                        ϕ_dist = Normal(ϕ,0.001));
chain = sample(test_model, NUTS(1_000,0.9), 1_000)
#check that STD of both vkick_div_vrel and log_m2 are small (sub 1% of median)
median_vkick_div_vrel = median(chain[:vkick_div_vrel])
std_vkick_div_vrel = std(chain[:vkick_div_vrel])
median_logm2_i = median(chain[:logm2_i])
std_logm2_i = std(chain[:logm2_i])
@test std_vkick_div_vrel < 0.01*median_vkick_div_vrel
@test std_logm2_i < 0.01*median_logm2_i
#verify that true solution is within the errors
@test abs(vkick_div_vrel-median_vkick_div_vrel) < std_vkick_div_vrel
@test abs(log10(m2_i)-median_logm2_i) < std_logm2_i

#verify that the same results are obtained if the order of the observables is shifted
test_model = SideKicks.createMCMCModel([:K1,:P,:K2,:e],
                         [K1,period,K2,eccentricity],
                         [K1,period,K2,eccentricity]*0.001,
                        logm1_i_dist = Normal(log10(m1_i),0.001),
                        frac_dist = truncated(Normal(frac,0.001),0.0001,0.9999),
                        cos_θ_dist = truncated(Normal(cos(θ),0.001),-1,1),
                        ϕ_dist = Normal(ϕ,0.001));
chain = sample(test_model, NUTS(1_000,0.9), 1_000)
#check that STD of both vkick_div_vrel and log_m2 are small (sub 1% of median)
median_vkick_div_vrel = median(chain[:vkick_div_vrel])
std_vkick_div_vrel = std(chain[:vkick_div_vrel])
median_logm2_i = median(chain[:logm2_i])
std_logm2_i = std(chain[:logm2_i])
@test std_vkick_div_vrel < 0.01*median_vkick_div_vrel
@test std_logm2_i < 0.01*median_logm2_i
#verify that true solution is within the errors
@test abs(vkick_div_vrel-median_vkick_div_vrel) < std_vkick_div_vrel
@test abs(log10(m2_i)-median_logm2_i) < std_logm2_i
