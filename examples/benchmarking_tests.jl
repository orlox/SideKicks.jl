using SideKicks
using BenchmarkTools


m_1i=10.0
m_2i=20.0
m_2f=15.0
vkick=10.0
θ=0.0
ϕ=0.0
a_i=1.0*au

##
@benchmark SideKicks.post_circular_kick_a(; m_1i=m_1i, m_2i=m_2i, m_2f=m_2f, vkick=vkick, θ=θ, ϕ=ϕ, a_i=a_i)

##

P_i = SideKicks.kepler_P_from_a(;m_1=m_1i, m_2=m_2i, a=a_i)
@benchmark SideKicks.post_circular_kick_P(; m_1i=m_1i, m_2i=m_2i, m_2f=m_2f, vkick=vkick, θ=θ, ϕ=ϕ, P_i=P_i)

##
@benchmark SideKicks.post_circular_kick_vsys(; m_1i=m_1i, m_2i=m_2i, m_2f=m_2f, vkick=vkick, θ=θ, ϕ=ϕ, a_i=a_i)

##
@benchmark SideKicks.post_circular_kick_parameters_timetest(; m_1i=m_1i, m_2i=m_2i, m_2f=m_2f, vkick=vkick, θ=θ, ϕ=ϕ, a_i=a_i)

##
@benchmark SideKicks.post_circular_kick_parameters_timetest_P(; m_1i=m_1i, m_2i=m_2i, m_2f=m_2f, vkick=vkick, θ=θ, ϕ=ϕ, P_i=SideKicks.kepler_P_from_a(m_1=m_1i, m_2=m_2i, a=a_i))
