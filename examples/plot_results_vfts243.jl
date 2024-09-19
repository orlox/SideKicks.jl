#=
# Plotting MCMC results from VFTS 243

In this example we produce corner plots from the results of VFTS 243 that were stored from
the previous example. We start by loading up 
=#

using CairoMakie
using SideKicks
using Distributions

results, observations, priors, metadata = SideKicks.ExtractResults("vfts243_results.hdf5", transpose_results=true)

##
#=
As an initial check, we can plot the quantities that were used as observations. A basic
consistency check is to verify these are consistent with the MCMC samples.
=#

plotting_props_obs_check = SideKicks.createPlottingProps([
    [:m1_f,    m_sun,    [15,40],        L"M_1\;[M_{\odot}]"],
    [:P_f,   day,      [10.35,10.45],    L"P_f\;[\mathrm{days}]"],
    [:e_f,   1,        [0,0.1],        L"e_f"],
    [:K1,    km_per_s, [77,90],        L"K_1  \;[\mathrm{km s}^{-1}]"],
])

cp = create_corner_plot(results, plotting_props_obs_check,
    supertitle="VFTS 243 - observables",
    dists_to_plot = Dict(
            :m1_f => Normal(25.0,2.3),
            :e_f => Normal(0.017,0.012),
            :K1 => Normal(81.4, 1.3),
            :P_f => Normal(10.4031, 0.01)
        )
    )
save("vfts243_observables.png", cp.fig)

cp.fig

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

println("Temporary issue with docs not showing derived quantities")

save("vfts243_derived.pdf", f)

f

##
#=
As you might see, the results are not very flattering, but this is a consequence of
the very small number of samples used.
=#
