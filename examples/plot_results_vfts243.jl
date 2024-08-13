#=
# Plotting MCMC results from VFTS 243

In this example we produce corner plots from the results of VFTS 243 that were stored from
the previous example. We start by loading up 
=#

using CairoMakie
using SideKicks

results, observations, priors, metadata = SideKicks.ExtractResults("vfts243_results.hdf5")

##
#=
As an initial check, we can plot the quantities that were used as observations. A basic
consistency check is to verify these are consistent with the MCMC samples.
=#

plotting_props_obs_check = SideKicks.createPlottingProps([
    [:m1,    m_sun,    [15,40],        L"M_1\;[M_{\odot}]"],
    [:P_f,   day,      [10.3,10.5],    L"P_f\;[\mathrm{days}]"],
    [:e_f,   1,        [0,0.1],        L"e_f"],
    [:K1,    km_per_s, [77,85],        L"K_1  \;[\mathrm{km s}^{-1}]"],
    [:vf_N,    km_per_s, [130,170],        L"v_N  \;[\mathrm{km s}^{-1}]"],
    [:vf_E,    km_per_s, [380,430],        L"v_E  \;[\mathrm{km s}^{-1}]"],
    [:vf_r,    km_per_s, [257,263],        L"v_r  \;[\mathrm{km s}^{-1}]"],
    [:ω_f,   degree, [0,360],        L"\omega_f  \;[\mathrm{rad}]"],
])

f = create_corner_plot(results, plotting_props_obs_check,
    tickfontsize=10 ,
    xticklabelrotation=pi/4, 
    show_CIs=true,
    rowcolgap=8,
    fraction_1D = 0.9,
    supertitle="VFTS 243 - observables",
    )
save("vfts243_observables.png", f)

f

##
#=
And, after verifying the samples do correspond to our observational constraints, we can analyze the
consequences for explosion itself.
=#

plotting_props = SideKicks.createPlottingProps([
    [:m2,     m_sun,    [0,25],         L"M_2  \;[M_{\odot}]"],
    [:dm2,    m_sun,    [0, 4],        L"ΔM_2  \;[M_{\odot}]"],
    [:P,      day,      [8,12],        L"P  \;[\mathrm{days}]"],
    [:vkick, km_per_s,  [0,50],         L"v_{kick}  \;[\mathrm{km s}^{-1}]"],
    [:vsys,  km_per_s, [0,50],        L"v_{\mathrm{sys}} \;[\mathrm{km s}^{-1}]"],
])

f = create_corner_plot(results, plotting_props,
    tickfontsize=10 ,
    xticklabelrotation=pi/4, 
    show_CIs=true,
    supertitle="VFTS 243 - derived quantities",
    fraction_1D = 0.9,
    )

save("vfts243_derived.png", f)

f

##
#=
As you might see, the results are not very flattering, but this is a consequence of
the very small number of samples used.
=#