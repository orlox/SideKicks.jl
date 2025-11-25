#=
# Plotting MCMC results from VFTS 243

In this example we produce corner plots from the results of VFTS 243 that were stored from
the previous example. We start by loading up 
=#

using CairoMakie
using SideKicks
using CornerPlotting
using Distributions

results, observations, priors, metadata = SideKicks.ExtractResults("vfts243_results.hdf5");

##
#=
As an initial check, we can plot the quantities that were used as observations. A basic
consistency check is to verify these are consistent with the MCMC samples.
=#

names = [:m1_f, :P_f, :e_f, :K1]
scaling = Dict( :m1_f => m_sun, 
                :P_f  => day,     
                :e_f  => 1,       
                :K1   => km_per_s
               )
labels  = Dict( :m1_f => L"M_{1f}\;[M_{\odot}]",
                :P_f  => L"P_f\;[\mathrm{days}]",
                :e_f  => L"e_f",
                :K1   => L"K_1  \;[\mathrm{km s}^{-1}]",
               )

set_theme!(CornerPlotting.default_theme())
cp = CornerPlotting.CornerPlot(results, names, scaling=scaling, labels=labels, nbins=10, nbins_contour=10)

CornerPlotting.plot_extra_1D_distribution(cp, :m1_f, Normal(25.0,2.3))
CornerPlotting.plot_extra_1D_distribution(cp, :e_f,  Normal(0.017,0.012))
CornerPlotting.plot_extra_1D_distribution(cp, :K1,   Normal(81.4, 1.3))
CornerPlotting.plot_extra_1D_distribution(cp, :P_f,  Normal(10.4031, 0.01))

save("vfts243_observables.png", cp.fig)

cp.fig

##
#=
And, after verifying the samples do correspond to our observational constraints, we can analyze the
consequences for explosion itself.
=#

names = [:m2_i, :dm2, :P_i, :vkick, :vsys]
scaling = Dict( :m2_i  => m_sun,    
                :dm2   => m_sun,    
                :P_i   => day,      
                :vkick => km_per_s, 
                :vsys  => km_per_s, 
               )
labels = Dict( :m2_i => L"M_{2f}\;[M_{\odot}]",
               :dm2   =>  L"ΔM_2  \;[M_{\odot}]",
               :P_i   =>  L"P_i \;[\mathrm{days}]",
               :vkick =>  L"v_{kick}  \;[\mathrm{km s}^{-1}]",
               :vsys  =>  L"v_{\mathrm{sys}} \;[\mathrm{km s}^{-1}]",
              )


set_theme!(CornerPlotting.default_theme())
cp = CornerPlotting.CornerPlot(results, names, scaling=scaling, labels=labels, nbins=10, nbins_contour=10)

save("vfts243_derived.png", cp.fig)

cp.fig

##
#=
As you might see, the results are not very flattering, but this is a consequence of
the very small number of samples used.
=#
