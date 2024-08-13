using SideKicks
using HDF5
using CairoMakie

system_id = last(split(@__DIR__, "/"))
results_file = String(@__DIR__) * "/results-medium.hdf5"

results, observations, priors, metadata = SideKicks.ExtractResults(results_file)

##

plotting_props_obs_check = SideKicks.createPlottingProps([
    [:m1,    m_sun,    [15,40],        L"M_1\;[M_{\odot}]"],
    [:P_f,   day,      [10.3,10.5],    L"P_f\;[\mathrm{days}]"],
    [:e_f,   1,        [0,0.1],        L"e_f"],
    [:K1,    km_per_s, [77,85],        L"K_1  \;[\mathrm{km s}^{-1}]"],
    [:vf_N,    km_per_s, [130,170],        L"v_N  \;[\mathrm{km s}^{-1}]"],
    #[:vf_E,    km_per_s, [380,430],        L"v_E  \;[\mathrm{km s}^{-1}]"],
    #[:vf_r,    km_per_s, [257,263],        L"v_r  \;[\mathrm{km s}^{-1}]"],
    #[:ω_f,   degree, [0,360],        L"\omega_f  \;[\mathrm{rad}]"],
])

    #[:vf_N,    km_per_s, [130,170],        L"v_N  \;[\mathrm{km s}^{-1}]"],
    #[:vf_E,    km_per_s, [380,430],        L"v_E  \;[\mathrm{km s}^{-1}]"],
    #[:vf_r,    km_per_s, [257,263],        L"v_r  \;[\mathrm{km s}^{-1}]"],

f = create_corner_plot(results, plotting_props_obs_check,
    observations=observations,
    tickfontsize=10 ,
    xticklabelrotation=3*pi/8, 
    show_CIs=true,
    nbins=50, # RTW - play with this, it may look better with larger samples
    rowcolgap=8,
    fraction_1D = 0.9,
    supertitle=system_id *" - observables",
    )
save(String(@__DIR__)*"/"*system_id*"_observables.png", f)

f

##

plotting_props = SideKicks.createPlottingProps([
    #[:m1,     m_sun,    [15,40],        L"M_1\;[M_{\odot}]"],
    #[:m2_f,   m_sun,    [6,14],         L"M_{2f}\;[M_{\odot}]"],
    #[:P_f,    day,      [10.37,10.43],  L"P_f\;[\mathrm{days}]"],
    #[:e_f,    1,        [0,0.1],        L"e_f"],
    #[:m2,     m_sun,    missing,         L"M_2  \;[M_{\odot}]"],
    [:m2,     m_sun,    [0,25],         L"M_2  \;[M_{\odot}]"],
    [:dm2,    m_sun,    [0, 4],        L"ΔM_2  \;[M_{\odot}]"],
    #[:dm2,    m_sun,    missing,        L"ΔM_2  \;[M_{\odot}]"],
    [:P,      day,      [8,12],        L"P  \;[\mathrm{days}]"],
    #[:a,      r_sun,    missing,        L"a  \;[R_{odot}]"],
    #[:i_f,    pi,       missing,        L"i_f  \;[π rad]"],
    [:vkick, km_per_s,  [0,50],         L"v_{kick}  \;[\mathrm{km s}^{-1}]"],
    #[:a_f,   r_sun,     missing,        L"a_f  \;[R_{odot}]"],
    #[:K1,    km_per_s,  missing,        L"K_1  \;[\mathrm{km s}^{-1}]"],
    #[:K2,    km_per_s,  missing,        L"K_2  \;[\mathrm{km s}^{-1}]"],
    #[:frac,  1,         missing,        L"f_{fb}"],
    [:vsys,  km_per_s, [0,50],        L"v_{\mathrm{sys}} \;[\mathrm{km s}^{-1}]"],
    #[:i_f,  degree, [0,90],        L"i_f"],
    #[:omega_f,  degree, [0,360],        L"\omega_f"],
    #[:Omega_f,  degree, [0,360],        L"\Omega_f"],
])

f = create_corner_plot(results, plotting_props,
    tickfontsize=10 ,
    xticklabelrotation=pi/4, 
    show_CIs=true,
    supertitle=system_id *" - derived quantities",
    fraction_1D = 0.9,
    )

save(String(@__DIR__)*"/"*system_id*"_derived.png", f)

f   

##

# Extra Eccentric plotting

plotting_props = SideKicks.createPlottingProps([
    #[:m1,     m_sun,    [15,40],        L"M_1\;[M_{\odot}]"],
    #[:m2_f,   m_sun,    [6,14],         L"M_{2f}\;[M_{\odot}]"],
    #[:P_f,    day,      [10.37,10.43],  L"P_f\;[\mathrm{days}]"],
    #[:e_f,    1,        [0,0.1],        L"e_f"],
    #[:m2,     m_sun,    missing,         L"M_2  \;[M_{\odot}]"],
    #[:m2,     m_sun,    [0,25],         L"M_2  \;[M_{\odot}]"],
    #[:dm2,    m_sun,    [0, 4],        L"ΔM_2  \;[M_{\odot}]"],
    #[:dm2,    m_sun,    missing,        L"ΔM_2  \;[M_{\odot}]"],
    #[:P,      day,      [8,12],        L"P  \;[\mathrm{days}]"],
    #[:a,      r_sun,    missing,        L"a  \;[R_{odot}]"],
    #[:i_f,    pi,       missing,        L"i_f  \;[π rad]"],
    [:vkick, km_per_s,  [0,50],         L"v_{kick}  \;[\mathrm{km s}^{-1}]"],
    #[:a_f,   r_sun,     missing,        L"a_f  \;[R_{odot}]"],
    #[:K1,    km_per_s,  missing,        L"K_1  \;[\mathrm{km s}^{-1}]"],
    #[:K2,    km_per_s,  missing,        L"K_2  \;[\mathrm{km s}^{-1}]"],
    #[:frac,  1,         missing,        L"f_{fb}"],
    [:vsys,  km_per_s, [0,50],        L"v_{\mathrm{sys}} \;[\mathrm{km s}^{-1}]"], 
    [:vf_N,  km_per_s, [0,50],        L"v_{\mathrm{N}} \;[\mathrm{km s}^{-1}]"], 
    [:vf_E,  km_per_s, [0,50],        L"v_{\mathrm{E}} \;[\mathrm{km s}^{-1}]"], 
    [:vf_r,  km_per_s, [0,50],        L"v_{\mathrm{r}} \;[\mathrm{km s}^{-1}]"], 
])

f = create_corner_plot(results, plotting_props,
    tickfontsize=10 ,
    xticklabelrotation=pi/4, 
    show_CIs=true,
    supertitle=system_id *" - derived quantities (ecc)"
    )

save(String(@__DIR__)*"/"*system_id*"_derived_ecc.png", f)

f   

##

# Plot everything

plotting_props = SideKicks.createPlottingProps([
    [:m1,     m_sun,    missing,  L"M_1\;[M_{\odot}]"],
    [:m2_f,   m_sun,    missing,  L"M_{2f}\;[M_{\odot}]"],
    [:P_f,    day,      missing,  L"P_f\;[\mathrm{days}]"],
    [:e_f,    1,        missing,  L"e_f"],
    [:m2,     m_sun,    missing,   L"M_2  \;[M_{\odot}]"],
    [:m2,     m_sun,    missing,  L"M_2  \;[M_{\odot}]"],
    [:dm2,    m_sun,    missing, L"ΔM_2  \;[M_{\odot}]"],
    [:dm2,    m_sun,    missing,  L"ΔM_2  \;[M_{\odot}]"],
    [:P,      day,      missing, L"P  \;[\mathrm{days}]"],
    [:a,      r_sun,    missing,  L"a  \;[R_{odot}]"],
    [:i_f,    pi,       missing,  L"i_f  \;[π rad]"],
    [:vkick, km_per_s,  missing,  L"v_{kick}  \;[\mathrm{km s}^{-1}]"],
    [:a_f,   r_sun,     missing,  L"a_f  \;[R_{odot}]"],
    [:K1,    km_per_s,  missing,  L"K_1  \;[\mathrm{km s}^{-1}]"],
    [:K2,    km_per_s,  missing,  L"K_2  \;[\mathrm{km s}^{-1}]"],
    [:frac,  1,         missing,  L"f_{fb}"],
    [:vsys,  km_per_s,  missing, L"v_{\mathrm{sys}} \;[\mathrm{km s}^{-1}]"], 
    [:vf_N,   km_per_s,  missing, L"v_{\mathrm{N}} \;[\mathrm{km s}^{-1}]"], 
    [:vf_E,   km_per_s,  missing, L"v_{\mathrm{E}} \;[\mathrm{km s}^{-1}]"], 
    [:vf_r,   km_per_s,  missing, L"v_{\mathrm{r}} \;[\mathrm{km s}^{-1}]"], 
])

f = create_corner_plot(results, plotting_props,
    tickfontsize=10 ,
    xticklabelrotation=pi/4, 
    show_CIs=true,
    supertitle=system_id *" - master plot"
    )

save(String(@__DIR__)*"/"*system_id*"_master.png", f)

f   

##
