using StatsBase
using CairoMakie
using CornerPlotting

export create_corner_plot 

"""
    mutable struct PlottingProps

PlottingProps contains the props, units, ranges, and names (latex encouraged)
for each property to plot.
"""

@kwdef mutable struct PlottingProps
    props::Vector{Symbol}
    units::Vector{Float64}
    ranges::Vector{Any}    # RTW TODO
    names::Vector{AbstractString}
end

function createPlottingProps(props_matrix::Vector{Vector{Any}})
    props_matrix = stack(props_matrix) # make into actual matrix
    return PlottingProps(
        props  = props_matrix[1,:],
        units  = props_matrix[2,:],
        ranges = props_matrix[3,:],
        names  = props_matrix[4,:])
end

function addPlottingProp(props_matrix::PlottingProps, new_prop::Vector{Any})
    props  = props_matrix.props
    units  = props_matrix.units 
    ranges = props_matrix.ranges 
    names  = props_matrix.names
    push!(props,  new_prop[1])
    push!(units , new_prop[2])
    push!(ranges, new_prop[3])
    push!(names,  new_prop[4])
    return PlottingProps(
        props  = props,
        units  = units,
        ranges = ranges,
        names  = names)
end

"""
    create_corner_plot(results, plotting_props; 
        observations=nothing, fig=Figure(), supertitle=nothing,
        fractions=[0.68,0.95,0.997], fraction_1D=0.68, 
        show_CIs=true, nbins=100, nbins_contour=20, rowcolgap=10, 
        xticklabelrotation=pi/4, labelfontsize=16, tickfontsize=10, supertitlefontsize=30)

Description
Function to create corner plot for selected (sub-)set of parameters from the MCMC output.

# Arguments:
- results:             the extracted results hdf5 object from a previous MCMC run                  
- plotting_props:      the plotting properties object containing which properties and ranges to plot                  
- observations:        any observations that should be included in the plots for comparison          
- fig:                 a figure, if needed
- supertitle:          the title of the plot          
- fractions:           the area fraction to determine different colored regions 
- fraction_1D:         the area fraction to include in the confidence interval bounds
- show_CIs:            whether to include confidence intervals
- nbins:               number of bins, identical for all parameters   
- nbins_contour:       number of bins for the contour plots
- rowcolgap:           spacing between the axes
- xticklabelrotation:  rotating (in rad) of the x-axis tick labels                
- labelfontsize:       fontsize of the parameter labels           
- tickfontsize:        fontsize of the tick labels
- supertitlefontsize:  fontsize of the title

# Output:
- fig:                 the newly created figure 
"""
function create_corner_plot(results, plotting_props; 
        dists_to_plot=nothing,
        fig=Figure(), supertitle=nothing,
        fraction_1D=0.9, fractions_2D=[0.9], 
        show_CIs=true, nbins=100, nbins_contour=30,
        supertitlefontsize=30, use_corner_plotting_theme=true)
 
    # TODO: Is there a way to add the priors? They are often modified
    # versions of the plotted parameters, this may be very non-trivial
    props::Vector{Symbol}  = [] 
    units  = Dict() 
    names  = Dict()
    ranges = Dict()
    # If extra plotting props are included that can't be used, scrap these and only plot the good ones. 
    available_props = keys(results)
    print_avail = false
    @show available_props
    for ii in eachindex(plotting_props.props)
        if plotting_props.props[ii] ∈ available_props
            prop = plotting_props.props[ii]
            push!(props,  prop)
            units[prop] = plotting_props.units[ii]
            names[prop] = plotting_props.names[ii]
            if !ismissing(plotting_props.ranges[ii])
                ranges[prop] = plotting_props.ranges[ii] 
            end
        else
            println("Prop "*string(plotting_props.props[ii])*" ignored")
            print_avail = true
        end
    end
    if print_avail
        println("Available props are "*string(available_props))
    end
    if use_corner_plotting_theme
        set_theme!(CornerPlotting.default_theme())
    end
    corner_plot = CornerPlotting.CornerPlot(results,props;labels=names,scaling=units)

    if !isnothing(dists_to_plot)
        for name in keys(dists_to_plot)
            CornerPlotting.plot_extra_1D_distribution(corner_plot, name, dists_to_plot[name])
        end
    end

    if !isnothing(supertitle)
        Label(corner_plot.fig[0,:], text=supertitle, fontsize=supertitlefontsize)
    end

    return corner_plot
end