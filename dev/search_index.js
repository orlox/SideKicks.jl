var documenterSearchIndex = {"docs":
[{"location":"2_plot_results/","page":"Plotting MCMC results from VFTS 243","title":"Plotting MCMC results from VFTS 243","text":"EditURL = \"../../examples/plot_results_vfts243.jl\"","category":"page"},{"location":"2_plot_results/#Plotting-MCMC-results-from-VFTS-243","page":"Plotting MCMC results from VFTS 243","title":"Plotting MCMC results from VFTS 243","text":"","category":"section"},{"location":"2_plot_results/","page":"Plotting MCMC results from VFTS 243","title":"Plotting MCMC results from VFTS 243","text":"In this example we produce corner plots from the results of VFTS 243 that were stored from the previous example. We start by loading up","category":"page"},{"location":"2_plot_results/","page":"Plotting MCMC results from VFTS 243","title":"Plotting MCMC results from VFTS 243","text":"using CairoMakie\nusing SideKicks\n\nresults, observations, priors, metadata = SideKicks.ExtractResults(\"vfts243_results.hdf5\")","category":"page"},{"location":"2_plot_results/","page":"Plotting MCMC results from VFTS 243","title":"Plotting MCMC results from VFTS 243","text":"As an initial check, we can plot the quantities that were used as observations. A basic consistency check is to verify these are consistent with the MCMC samples.","category":"page"},{"location":"2_plot_results/","page":"Plotting MCMC results from VFTS 243","title":"Plotting MCMC results from VFTS 243","text":"plotting_props_obs_check = SideKicks.createPlottingProps([\n    [:m1_f,    m_sun,    [15,40],        L\"M_1\\;[M_{\\odot}]\"],\n    [:P_f,   day,      [10.35,10.45],    L\"P_f\\;[\\mathrm{days}]\"],\n    [:e_f,   1,        [0,0.1],        L\"e_f\"],\n    [:K1,    km_per_s, [77,90],        L\"K_1  \\;[\\mathrm{km s}^{-1}]\"],\n])\n\nf = create_corner_plot(results, plotting_props_obs_check,\n    supertitle=\"VFTS 243 - observables\",\n    )\nsave(\"vfts243_observables.png\", f)\n\nf","category":"page"},{"location":"2_plot_results/","page":"Plotting MCMC results from VFTS 243","title":"Plotting MCMC results from VFTS 243","text":"And, after verifying the samples do correspond to our observational constraints, we can analyze the consequences for explosion itself.","category":"page"},{"location":"2_plot_results/","page":"Plotting MCMC results from VFTS 243","title":"Plotting MCMC results from VFTS 243","text":"plotting_props = SideKicks.createPlottingProps([\n    [:m2,     m_sun,    [0,25],         L\"M_2  \\;[M_{\\odot}]\"],\n    [:dm2,    m_sun,    [0, 4],        L\"ΔM_2  \\;[M_{\\odot}]\"],\n    [:P,      day,      [8,12],        L\"P  \\;[\\mathrm{days}]\"],\n    [:vkick, km_per_s,  [0,50],         L\"v_{kick}  \\;[\\mathrm{km s}^{-1}]\"],\n    [:vsys,  km_per_s, [0,50],        L\"v_{\\mathrm{sys}} \\;[\\mathrm{km s}^{-1}]\"],\n])\n\nf = create_corner_plot(results, plotting_props,\n    tickfontsize=10 ,\n    xticklabelrotation=pi/4,\n    show_CIs=true,\n    supertitle=\"VFTS 243 - derived quantities\",\n    fraction_1D = 0.9,\n    )\n\nprintln(\"Temporary issue with docs not showing derived quantities\")\n\n#save(\"vfts243_derived.pdf\", f)\n\n#f","category":"page"},{"location":"2_plot_results/","page":"Plotting MCMC results from VFTS 243","title":"Plotting MCMC results from VFTS 243","text":"As you might see, the results are not very flattering, but this is a consequence of the very small number of samples used.","category":"page"},{"location":"2_plot_results/","page":"Plotting MCMC results from VFTS 243","title":"Plotting MCMC results from VFTS 243","text":"","category":"page"},{"location":"2_plot_results/","page":"Plotting MCMC results from VFTS 243","title":"Plotting MCMC results from VFTS 243","text":"This page was generated using Literate.jl.","category":"page"},{"location":"1_run_inference/","page":"Inference example for VFTS 243","title":"Inference example for VFTS 243","text":"EditURL = \"../../examples/run_inference_vfts243.jl\"","category":"page"},{"location":"1_run_inference/#Inference-example-for-VFTS-243","page":"Inference example for VFTS 243","title":"Inference example for VFTS 243","text":"","category":"section"},{"location":"1_run_inference/","page":"Inference example for VFTS 243","title":"Inference example for VFTS 243","text":"In this example we try to match the observed properties of VFTS 243. We start by loading up the SideKicks package, as well the Distributions package which we use to define the priors.","category":"page"},{"location":"1_run_inference/","page":"Inference example for VFTS 243","title":"Inference example for VFTS 243","text":"using SideKicks\nusing Distributions","category":"page"},{"location":"1_run_inference/","page":"Inference example for VFTS 243","title":"Inference example for VFTS 243","text":"We set up our observations using the @Observations macro. This macro will return a tuple containing an instance of Observation, as well as a string with the code needed to initialize the Observation instance. The reason for this is that when we save results, we also want to save the precise observations that were used for the MCMC.","category":"page"},{"location":"1_run_inference/","page":"Inference example for VFTS 243","title":"Inference example for VFTS 243","text":"Currently, The likelihoods are taken to be Normal distributions. This can potentially be generalized in the future to arbitrary user-defined likelihood functions.","category":"page"},{"location":"1_run_inference/","page":"Inference example for VFTS 243","title":"Inference example for VFTS 243","text":"obs = SideKicks.@Observations([\n    [:P_f,  10.4031, 0.01,   day],\n    [:e_f,  0.017,   0.012,  1],\n    [:m1_f, 25.0,    2.3,    m_sun],\n    [:K1,   81.4,    1.3,    km_per_s],\n    [:v_N,  138.8,    7.6,     km_per_s], # Gaia 4' w/ MCMC\n    [:v_E,  409.3,    9.3,     km_per_s], # Gaia 4' w/ MCMC\n    [:v_r,  261.5,   0.42,    km_per_s], # Almeida\n    [:ω_f,  66,      53,     degree]\n])","category":"page"},{"location":"1_run_inference/","page":"Inference example for VFTS 243","title":"Inference example for VFTS 243","text":"Next, we choose our priors. In here we are using distributions defined in the Distributions.jl` package.","category":"page"},{"location":"1_run_inference/","page":"Inference example for VFTS 243","title":"Inference example for VFTS 243","text":"priors = SideKicks.@Priors(\n    logm1_dist = Uniform(0.1,3), # in log(Msun)\n    logm2_dist = Uniform(0.1,3), # in log(Msun)\n    logP_dist  = Uniform(-1,3),   # in log(days)\n    vkick_dist = Uniform(0,4), # in 100 km/s\n    frac_dist  = Uniform(0,1.0),\n    e_dist = Uniform(0,0.01),\n    venv_N_100kms_dist = Normal(146/100, 40/100), # Gaia 4' w/ MCMC\n    venv_E_100kms_dist = Normal(393/100, 42/100), # Gaia 4' w/ MCMC\n    venv_r_100kms_dist = Normal(270.3/100, 11.1/100) # Almeida w/ MCMC\n)","category":"page"},{"location":"1_run_inference/","page":"Inference example for VFTS 243","title":"Inference example for VFTS 243","text":"Finally, we run the MCMC using the NUTS algorithm. In the example below we are computing 3 chains with 100 samples each. This in practice is very low, and a ballpark suggestion would be to total at least sim10^4 samples for any meaningful analysis, and an order of magnitude higher for final science runs.","category":"page"},{"location":"1_run_inference/","page":"Inference example for VFTS 243","title":"Inference example for VFTS 243","text":"kick_mcmc = SideKicks.KickMCMC(\n        which_model = :general,\n        observations = obs,\n        priors = priors,\n        nuts_warmup_count = 200,\n        nuts_acceptance_rate = 0.8,\n        nsamples = 200,\n        nchains = 4)","category":"page"},{"location":"1_run_inference/","page":"Inference example for VFTS 243","title":"Inference example for VFTS 243","text":"Results from the MCMC can be saved to an HDF5 file. In practice, it is ideal to separate the computation of the MCMC and its analysis, so plotting of the results is done separetely.","category":"page"},{"location":"1_run_inference/","page":"Inference example for VFTS 243","title":"Inference example for VFTS 243","text":"SideKicks.SaveResults(\"vfts243_results.hdf5\", kick_mcmc)","category":"page"},{"location":"1_run_inference/","page":"Inference example for VFTS 243","title":"Inference example for VFTS 243","text":"","category":"page"},{"location":"1_run_inference/","page":"Inference example for VFTS 243","title":"Inference example for VFTS 243","text":"This page was generated using Literate.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = SideKicks","category":"page"},{"location":"#SideKicks","page":"Home","title":"SideKicks","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for SideKicks, a Julia package for analyzing the observed properties of stellar binaries containing a compact object, in order to perform parameter inference on the supernova mass loss and natal kick. ","category":"page"},{"location":"#Package-features","page":"Home","title":"Package features","text":"","category":"section"},{"location":"##-TODO-RTW","page":"Home","title":"# TODO - RTW","text":"","category":"section"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [SideKicks]","category":"page"},{"location":"#SideKicks.SideKicks","page":"Home","title":"SideKicks.SideKicks","text":"Main module for SideKicks.jl – a analysis package for performing parameter inference on compact object stellar binaries.\n\n\n\n\n\n","category":"module"},{"location":"#SideKicks.KickMCMC","page":"Home","title":"SideKicks.KickMCMC","text":"mutable struct KickMCMC\n\nKickMCMC contains the MCMC model and Observations structs, the Results dict, the chains that resulted from the MCMC, and the parameters that went into the sampler.\n\n\n\n\n\n","category":"type"},{"location":"#SideKicks.ModVonMises","page":"Home","title":"SideKicks.ModVonMises","text":"struct ModVonMises{T1<:Real, T2<:Real} <: ContinuousUnivariateDistribution\n\nThis is just a wrapper on top of the VonMises distribution (as defined in Distributions.jl) to extend its domain. This is because the domain of VonMises is defined to be [μ-π, μ+π], and the angles we are concerned with range from [0,2π]\n\n\n\n\n\n","category":"type"},{"location":"#SideKicks.Observations","page":"Home","title":"SideKicks.Observations","text":"mutable struct Observations\n\nObservations contains the symbols, values, errors, and units of each observed parameter.\n\n\n\n\n\n","category":"type"},{"location":"#SideKicks.Priors","page":"Home","title":"SideKicks.Priors","text":"mutable struct Priors\n\nPriors contains the prior distribution of each of the desired parameters\n\n\n\n\n\n","category":"type"},{"location":"#SideKicks.WrappedCauchy","page":"Home","title":"SideKicks.WrappedCauchy","text":"struct WrappedCauchy{T1<:Real, T2<:Real} <: ContinuousUnivariateDistribution\n\nThe WrappedCauchy distribution resembles the Cauchy distribution defined on the unit  circle from 0 to 2π, with the endpoints wrapped back to each other.\n\n\n\n\n\n","category":"type"},{"location":"#SideKicks.RV_semiamplitude_K1-Tuple{}","page":"Home","title":"SideKicks.RV_semiamplitude_K1","text":"RV_semiamplitude_K1(;m1, m2, P, e, i)\n\nCompute the amplitude of radial velocity variations given orbital parameters and masses\n\nArguments:\n\nm1:   mass of observed star             [g]\nm2:   mass of companion                 [g]\nP:    orbital period                    [s]\ne:    orbital eccentricity              [-]\ni:    orbital inclination               [rad]\n\nOutput:\n\nK1: amplitude of radial velocity variation of star 1  [cm/s]\n\n\n\n\n\n","category":"method"},{"location":"#SideKicks.arbitraryEjectaBH-Tuple{Any, Any}","page":"Home","title":"SideKicks.arbitraryEjectaBH","text":"arbitraryEjectaBH(m2_i, frac)\n\n#TODO Description\n\nArguments:\n\n#TODO\n\nm2_i:\nfrac:\n\nOutput:\n\n#TODO\n\n\n\n\n\n\n\n","category":"method"},{"location":"#SideKicks.create_1D_density-NTuple{6, Any}","page":"Home","title":"SideKicks.create_1D_density","text":"create_1D_density(axis, values, range, chain_weights, fraction_1D, nbins; color, linewidth)\n\nDescription Make the 1D density plots given the parameter values, ranges, and weights.\n\nArguments:\n\naxis:           the axis to make the plot\nvalues:         the values for the x-coordinate        \nrange:          the ranges for the x-coordinate        \nchain_weights:  the sample weighting from the MCMC\nfraction_1D:    the fractional area from which to compute the confidence intervals\nnbins:          number of bins, identical for all parameters   \ncolor:          the color of the density curve\nlinewidth:      the linewidth of the density curve\n\nOutput:\n\nx:              the x-coordinates of the density plot\nh:              the heights of the density plot\ny:              the normalized heights of the density plot\n\n\n\n\n\n","category":"method"},{"location":"#SideKicks.create_2D_density-NTuple{7, Any}","page":"Home","title":"SideKicks.create_2D_density","text":"create_2D_density(axis, values1, ranges1, values2, ranges2, chain_weights, fractions, nbins)\n\nDescription Make the 2D density plots given the parameter values, ranges, and weights.\n\nTODO: check that x- and y- descripters are correct, here and in below functions\n\nArguments:\n\naxis:           the axis to make the plot\nvalues1:        the values for the x-coordinate        \nranges1:        the ranges for the x-coordinate        \nvalues2:        the values for the y-coordinate        \nranges2:        the ranges for the y-coordinate        \nchain_weights:  the sample weighting from the MCMC\nfractions:      area fractions for defining contours\nnbins:          number of bins, identical for all parameters   \n\n\n\n\n\n","category":"method"},{"location":"#SideKicks.create_compound_1D_densities-NTuple{6, Any}","page":"Home","title":"SideKicks.create_compound_1D_densities","text":"create_compound_1D_densities(axis, values_matrix, range, chain_weights_matrix, fraction_1D, nbins)\n\n#TODO Description\n\nArguments:\n\naxis:                 the axis to make the plot\nvalues_matrix:        the values for each chain of the parameter \nrange:                the ranges for the x-coordinate        \nchainweightsmatrix: the sample weighting for each chain \nfraction_1D:          the fractional area from which to compute the confidence intervals\nnbins:                number of bins, identical for all parameters   \n\nOutput:\n\nxmin:                 the left boundary of the fraction_1D area interval\nxmode:                the mode of the data, within the provided range\nxmax:                 the right boundary of the fraction_1D area interval\n\n\n\n\n\n","category":"method"},{"location":"#SideKicks.create_corner_plot-Tuple{Any, Any}","page":"Home","title":"SideKicks.create_corner_plot","text":"create_corner_plot(results, plotting_props; \n    observations=nothing, fig=Figure(), supertitle=nothing,\n    fractions=[0.68,0.95,0.997], fraction_1D=0.68, \n    show_CIs=true, nbins=100, nbins_contour=20, rowcolgap=10, \n    xticklabelrotation=pi/4, labelfontsize=16, tickfontsize=10, supertitlefontsize=30)\n\nDescription Function to create corner plot for selected (sub-)set of parameters from the MCMC output.\n\nArguments:\n\nresults:             the extracted results hdf5 object from a previous MCMC run                  \nplotting_props:      the plotting properties object containing which properties and ranges to plot                  \nobservations:        any observations that should be included in the plots for comparison          \nfig:                 a figure, if needed\nsupertitle:          the title of the plot          \nfractions:           the area fraction to determine different colored regions \nfraction_1D:         the area fraction to include in the confidence interval bounds\nshow_CIs:            whether to include confidence intervals\nnbins:               number of bins, identical for all parameters   \nnbins_contour:       number of bins for the contour plots\nrowcolgap:           spacing between the axes\nxticklabelrotation:  rotating (in rad) of the x-axis tick labels                \nlabelfontsize:       fontsize of the parameter labels           \ntickfontsize:        fontsize of the tick labels\nsupertitlefontsize:  fontsize of the title\n\nOutput:\n\nfig:                 the newly created figure \n\n\n\n\n\n","category":"method"},{"location":"#SideKicks.create_general_mcmc_model-Tuple{}","page":"Home","title":"SideKicks.create_general_mcmc_model","text":"create_general_mcmc_model(observations, observed_values, observed_errors)\n\nCreate a Turing model to perform an MCMC sampling of the pre-explosion  and kick properties of a system, assuming pre-explosion eccentricity.\n\nRTW does acos just work? Do I need to worry about domain/range issues?\n\nCheck velocities, the conversions are a bit funky\n\nArguments:\n\nobservations:    the parameters taken from observations [Vector{Symbol}]\nobserved_values: the values of the parameters           [Vector{Float64}] \nobserved_errors: the errors of the observations         [Vector{Float64}]\n\nOutput:\n\nkickmodel: A Turing model for sampling\n\n\n\n\n\n","category":"method"},{"location":"#SideKicks.create_simplified_mcmc_model-Tuple{}","page":"Home","title":"SideKicks.create_simplified_mcmc_model","text":"create_simplified_mcmc_model(observations, observed_values, observed_errors)\n\nDescription Create a Turing model to perform an MCMC sampling of the pre-explosion  and kick properties of a system, assuming pre-explosion circularity.\n\nRTW this is more simplistic than just using a circular model, it's also assuming you know the eccentricity and don't care about radial velocity etc.\n\nArguments:\n\nobservations:    the parameters taken from observations [Vector{Symbol}]\nobserved_values: the values of the parameters           [Vector{Float64}] \nobserved_errors: the errors of the observations         [Vector{Float64}]\n\nOutput:\n\nkickmodel: A Turing model for sampling\n\n\n\n\n\n","category":"method"},{"location":"#SideKicks.get_bounds_for_fractions-Tuple{Any, Any}","page":"Home","title":"SideKicks.get_bounds_for_fractions","text":"get_bounds_for_fractions(h, fractions)\n\nDescription Calculate the bounds containing the specified fraction(s) of area.\n\nArguments:\n\nh:         the densities contained in the bins\nfractions: the fractional area that should be bounded\n\nOutput:\n\nbounds:    the limits of the bounding area\n\n\n\n\n\n","category":"method"},{"location":"#SideKicks.kepler_P_from_a-Tuple{}","page":"Home","title":"SideKicks.kepler_P_from_a","text":"kepler_P_from_a(;m1, m2, a)\n\nObtain period from semimajor axis using Kepler's third law\n\nArguments:\n\nm1: mass of first companion      [g]\nm2: mass of 2nd companion        [g]\na:  semi-major axis of the orbit [cm]\n\nOutput:\n\nP: the orbital period            [s]\n\n\n\n\n\n","category":"method"},{"location":"#SideKicks.kepler_a_from_P-Tuple{}","page":"Home","title":"SideKicks.kepler_a_from_P","text":"kepler_a_from_P(;m1, m2, P)\n\nObtain semimajor axis from period using Kepler's third law\n\nArguments:\n\nm1: mass of first companion [g]\nm2: mass of 2nd companion   [g]\nP:  orbital period          [s]\n\nOutput:\n\na: semi-major axis of the orbit [cm]\n\n\n\n\n\n","category":"method"},{"location":"#SideKicks.post_supernova_circular_orbit_P-Tuple{}","page":"Home","title":"SideKicks.post_supernova_circular_orbit_P","text":"post_supernova_circular_orbit_P(;m1_i, m2_i, P_i, m1_f=-1, m2_f, vkick=0, θ=0, ϕ=0, vimp=0)\n\nSame as post_supernova_circular_orbit_a, except that it receives the initial orbital period as input and returns the final orbital period and eccentricity.\n\nArguments:\n\nm1_i:  pre-explosion  mass of non-exploding component [g]           \nm2_i:  pre-explosion  mass of exploding component     [g]       \nP_i:   pre-explosion orbital period                   [d]\nm1_f:  post-explosion mass of non-exploding component [g]           \nm2_f:  post-explosion mass of exploding component     [g]   \nvkick: kick velocity                                  [cm/s] \nθ:     polar kick angle (away from e_par)             [rad]\nϕ:     azimuthal kick angle (off of e_perp)           [rad]\nvimp:  imparted kick velocity on companion            [cm/s]     \n\nOutput:\n\nP_f: post-explosion orbital period                    [d]\ne_f: post-explosion excentricity                      [-]\n\n\n\n\n\n","category":"method"},{"location":"#SideKicks.post_supernova_circular_orbit_a-Tuple{}","page":"Home","title":"SideKicks.post_supernova_circular_orbit_a","text":"post_supernova_circular_orbit_a(;m1_i, m2_i, a_i, m1_f=-1.0, m2_f, vkick=0.0, θ=0.0, ϕ=0.0, vimp=0.0)\n\nCompute post-kick properties for a circular pre-explosion orbit. Equivalent to Tauris et al. (1999): Monthly Notices of the Royal Astronomical Society, Volume 310, Issue 4, pp. 1165-1169.\n\nArguments:\n\nm1_i:  pre-explosion  mass of non-exploding component [g]           \nm2_i:  pre-explosion  mass of exploding component     [g]       \na_i:   pre-explosion orbital separation               [cm]\nm1_f:  post-explosion mass of non-exploding component [g]           \nm2_f:  post-explosion mass of exploding component     [g]   \nvkick: kick velocity                                  [cm/s] \nθ:     polar kick angle (away from e_par)             [rad]\nϕ:     azimuthal kick angle (off of e_perp)           [rad]\nvimp:  imparted kick velocity on companion            [cm/s]     \n\nOutput:\n\na_f: post-explosion orbital separation                [cm]\ne_f: post-explosion excentricity                      [-]\n\n\n\n\n\n","category":"method"},{"location":"#SideKicks.post_supernova_circular_orbit_vsys-Tuple{}","page":"Home","title":"SideKicks.post_supernova_circular_orbit_vsys","text":"post_supernova_circular_orbit_vsys(;m1_i, m2_i, a_i, m1_f=-1, m2_f, vkick=0, θ=0, ϕ=0, vimp=0)\n\nCompute post-kick systemic velocity for a circular orbit Tauris et al. (1999): Monthly Notices of the Royal Astronomical Society, Volume 310, Issue 4, pp. 1165-1169.\n\nArguments:\n\nm1_i:  pre-explosion  mass of non-exploding component   [g]           \nm2_i:  pre-explosion  mass of exploding component       [g]       \na_i:   pre-explosion orbital separation                 [cm]\nm1_f:  post-explosion mass of non-exploding component [g]           \nm2_f:  post-explosion mass of exploding component     [g]   \nvkick: kick velocity                                  [cm/s] \nθ:     polar kick angle (away from e_par)             [rad]\nϕ:     azimuthal kick angle (off of e_perp)           [rad]\nvimp:  imparted kick velocity on companion            [cm/s]     \n\nOutput:\n\nvsys_f: post-explosion systemic velocity              [cm/s]\n\n\n\n\n\n","category":"method"},{"location":"#SideKicks.post_supernova_general_orbit_parameters-Tuple{}","page":"Home","title":"SideKicks.post_supernova_general_orbit_parameters","text":"post_supernova_general_orbit_parameters(;m1_i, m2_i, a_i, e_i=0, m1_f=-1, m2_f, vkick=0, θ=0, ϕ=0, vimp=0,\n    ν_i=0, Ω_i=0, ω_i=0, i_i=0)\n\nCompute post-kick properties for a general pre-explosion orbit  using equations from [Marchant, Willcox, Vigna-Gomez] TODO\n\nArguments:\n\nm1_i:  pre-explosion  mass of non-exploding component  [g]           \nm2_i:  pre-explosion  mass of exploding component      [g]       \na_i:   pre-explosion orbital separation                [cm]\ne_i:   pre-explosion orbital eccentricity              [-]\nm1_f:  post-explosion mass of non-exploding component  [g]           \nm2_f:  post-explosion mass of exploding component      [g]   \n\nvkick: kick velocity                                   [cm/s] \nθ:     polar kick angle (away from e_par)              [rad]\nϕ:     azimuthal kick angle (off of e_perp)            [rad]\nvimp:  imparted kick velocity on companion             [cm/s]     \nInitial orbital orientation angles: \nν_i: true anomaly                                    [rad]\nΩ_i: pre-explosion longitude of the ascending node   [rad]\nω_i: pre-explosion argument of periastron            [rad]\ni_i: pre-explosion inclination                       [rad]\n\nOutput: RTW: check!\n\na_f:   post-explosion orbital separation               [cm]\ne_f:   post-explosion orbital eccentricity             [-]\nΩ_f:   post-explosion longitude of ascending node      [rad]      \nω_f:   post-explosion argument of periastron           [rad]    \ni_f:   post-explosion inclination                      [rad]     \nv_n:   post-explosion systemic velocity, toward N      [rad]\nv_w:   post-explosion systemic velocity, toward W      [rad]      \nv_rad: post-explosion radial velocity, toward O        [rad]      \n\n\n\n\n\n","category":"method"},{"location":"#SideKicks.relative_velocity-Tuple{}","page":"Home","title":"SideKicks.relative_velocity","text":"relative_velocity(;m1, m2, a)\n\nCalculate the relative orbital velocity for a circular orbit.\n\nArguments:\n\nm1: mass of first companion        [g]\nm2: mass of 2nd companion          [g]\na:  semi-major axis of the orbit   [cm]\n\nOutput:\n\nv_rel: the relative velocity [cm/s]\n\n\n\n\n\n","category":"method"},{"location":"#SideKicks.restrictedEjectaBH-Tuple{Any, Any}","page":"Home","title":"SideKicks.restrictedEjectaBH","text":"restrictedEjectaBH(m2_i, frac; Ma=10, Mb=15, max_frac=1.0, min_frac=0.1)\n\n#TODO Description\n\nArguments:\n\n#TODO\n\nm2_i:\nfrac:\nMa=10:\nMb=15:\nmax_frac=1.0:\nmin_frac=0.1:\n\nOutput:\n\n#TODO\n\n\n\n\n\n\n\n","category":"method"}]
}
