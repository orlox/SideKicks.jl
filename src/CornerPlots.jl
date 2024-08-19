using StatsBase
using CairoMakie

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
        show_CIs=true, nbins=100, rowcolgap=10, 
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
- rowcolgap:           spacing between the axes
- xticklabelrotation:  rotating (in rad) of the x-axis tick labels                
- labelfontsize:       fontsize of the parameter labels           
- tickfontsize:        fontsize of the tick labels
- supertitlefontsize:  fontsize of the title

# Output:
- fig:                 the newly created figure 
"""

function create_corner_plot(results, plotting_props; 
        observations=nothing,
        fig=Figure(), supertitle=nothing,
        fractions_2D=[0.68,0.95,0.997], fraction_1D=0.68, 
        show_CIs=true, nbins=100,
        rowcolgap=10, xticklabelrotation=pi/4,
        labelfontsize=16, tickfontsize=10, supertitlefontsize=30)
 
    # TODO: Is there a way to add the priors? They are often modified
    # versions of the plotted parameters, this may be very non-trivial
    props  = [] 
    units  = [] 
    names  = [] 
    ranges = [] 
    # If extra plotting props are included that can't be used, scrap these and only plot the good ones. 
    available_props = keys(results)
    print_avail = false
    for ii in eachindex(plotting_props.props)
        if plotting_props.props[ii] ∈ available_props
            push!(props,  plotting_props.props[ii])
            push!(units,  plotting_props.units[ii])
            push!(names,  plotting_props.names[ii])
            push!(ranges, plotting_props.ranges[ii]) 
        else
            println("Prop "*string(plotting_props.props[ii])*" ignored")
            print_avail = true
        end
    end
    if print_avail
        println("Available props are "*string(available_props))
    end
    # Add ranges if none supplied
    num_props = length(props)
    for ii in 1:num_props
        values = vec(results[props[ii]])/units[ii]
        minval = minimum(values)
        maxval = maximum(values)
        if ismissing(ranges[ii]) || minval > ranges[ii][2] || maxval < ranges[ii][1]
            ranges[ii] = [minval, maxval]
            println("Range set for ", props[ii], ": ", ranges[ii])
        end
    end
    # Create 2D density plots
    num_col = num_props
    for ii in 1:num_col-1
        for  jj in ii+1:num_col
            axis = Axis(fig[jj+1,ii], xtickalign=1, xtickcolor = :black, ytickalign=1, ytickcolor = :black, 
                        aspect=1,
                        xlabel=names[ii], ylabel=names[jj], 
                        xlabelsize=labelfontsize, ylabelsize=labelfontsize,
                        xticklabelrotation=xticklabelrotation,
                        xticklabelsize=tickfontsize, yticklabelsize=tickfontsize)
            create_2D_density(axis, vec(results[props[ii]])/units[ii], ranges[ii], vec(results[props[jj]])/units[jj], ranges[jj], vec(results[:weights]), fractions_2D, nbins)
            #create_2D_density(axis, vec(results[props[ii]])/units[ii], ranges[ii], vec(results[props[jj]])/units[jj], ranges[jj], vec(results[:weights]), nbins)
            if ii>1
                hideydecorations!(axis, ticks=false, minorticks=false)
            end
            if jj!=num_col
                hidexdecorations!(axis,ticks=false, minorticks=false)
            end         
        end  
    end 
    # Create 1D PDFs along the diagonal
    latex_bounds_array = Array{AbstractString}(undef, num_col)
    for ii in 1:num_col
        axis = Axis(fig[ii+1,ii], xgridvisible = false, ygridvisible = false, xtickalign=1, xlabel=names[ii],
                    aspect=1,
                        xlabelsize=labelfontsize, ylabelsize=labelfontsize,
                        xticklabelrotation=xticklabelrotation,
                        xticklabelsize=tickfontsize, yticklabelsize=tickfontsize)
        (xmin, xmode, xmax) = create_compound_1D_densities(axis, results[props[ii]]/units[ii], ranges[ii], results[:weights], fraction_1D, nbins)
        
        # Add observations
        if !isnothing(observations)
            if props[ii] ∈ observations.props
                println("adding obs for "*String(props[ii]))
                # find index of prop and get mean and std
                idx = findall(x->x==props[ii], observations.props)[1]
                mean = observations.vals[idx]
                errs = observations.errs[idx]
                xarr = LinRange(ranges[ii][1], ranges[ii][2], 100)
                lines!(axis, xarr, pdf(Normal(mean, errs), xarr), color=(:red, 0.4), linewidth=3) # linestyle=:dot, 
            end
        end

        # Remove labels on the diagonals, except xlabel on the bottom right
        hideydecorations!(axis)
        if ii !=num_col
            hidexdecorations!(axis,ticks=false, minorticks=false)
        end
        # Configure confidence intervals
        str_xmode = round(xmode, sigdigits=3)
        str_upp = round(xmax-xmode, sigdigits=3)
        str_low = round(xmode-xmin, sigdigits=3)
        latex_bounds = L"%$(str_xmode)^{+%$(str_upp)}_{-%$(str_low)}"
        latex_bounds_array[ii] = latex_bounds
        println(names[ii]*"="*latex_bounds)
        if show_CIs
            if ii == 1
                Label(fig[ii,ii], latex_bounds, 
                      valign=:bottom, fontsize=20,
                      tellwidth=false)
                aspect=3
            else
                Label(fig[ii,ii], latex_bounds, 
                      valign=:bottom, fontsize=20,
                      tellwidth=false, tellheight=false)
                aspect=1
            end
            axis = Axis(fig[ii,ii], aspect=aspect,
                xgridvisible = false, ygridvisible = false)
            hidedecorations!(axis)
            hidespines!(axis)
        end
    end     

    # Readjust rows and columns
    rowgap!(fig.layout, rowcolgap)
    colgap!(fig.layout, rowcolgap)

    if !isnothing(supertitle)
        Label(fig[0,:], text=supertitle, fontsize=supertitlefontsize)
    end

    resolution=600 .*(1,1)
    resize!(fig.scene, resolution)
    return fig
end

"""
    create_2D_density(axis, values1, ranges1, values2, ranges2, chain_weights, fractions, nbins)

Description
Make the 2D density plots given the parameter values, ranges, and weights.

# TODO: check that x- and y- descripters are correct, here and in below functions

# Arguments:
- axis:           the axis to make the plot
- values1:        the values for the x-coordinate        
- ranges1:        the ranges for the x-coordinate        
- values2:        the values for the y-coordinate        
- ranges2:        the ranges for the y-coordinate        
- chain_weights:  the sample weighting from the MCMC
- fractions:      area fractions for defining contours
- nbins:          number of bins, identical for all parameters   
"""
function create_2D_density(axis, values1, ranges1, values2, ranges2, chain_weights, fractions, nbins)

    filter1 = values1 .> ranges1[1] .&& values1 .< ranges1[2] 
    filter2 = values2 .> ranges2[1] .&& values2 .< ranges2[2]
    # RTW: make this check more robust
    if sum(filter1) == 0
        println("problem with 1")
        println(ranges1, " ", minimum(values1), " ", maximum(values1))
        return
    end
    if sum(filter2) == 0
        println("problem with 2")
        println(ranges2, " ", minimum(values2), " ", maximum(values2))
        return
    end
    filter = filter1 .&& filter2
    values1 = values1[filter]
    values2 = values2[filter]
    chain_weights = weights(chain_weights[filter]) # weights is a StatsBase function
    
    h = fit(Histogram, (values1, values2), chain_weights, nbins=nbins)
    x = (h.edges[2][2:end] .+ h.edges[2][1:end-1])./2
    y = (h.edges[1][2:end] .+ h.edges[1][1:end-1])./2
    bounds = get_bounds_for_fractions(h, fractions)
    contour_matrix = zeros(size(h.weights))
    nbounds = length(bounds)
    for ii in 1:nbounds
        contour_matrix[h.weights .> bounds[ii]] .= (ii-1)/(nbounds-1)
    end
    #colormap=(Reverse(:dense), 0.2) # 2nd number is alpha
    colormap=(:dense, 1.0) # 2nd number is alpha
    #contourf!(axis, y, x, contour_matrix, colormap=colormap, colorrange=(0,1), extendlow=:white)
    println(maximum(h.weights))
    println(minimum(h.weights))
    println(bounds)
    #contourf!(axis, h.weights, levels=bounds, colormap=colormap, 
    #heatmap!(axis, y, x, h.weights, colormap=:dense)
    heatmap!(axis, y, x, h.weights, colormap=colormap)
    colormap = [:red, :red, :red]
    contour!(axis, y, x, h.weights, levels=bounds, colormap=colormap, linewidth=1)
              #extendlow=:white)#, extendhigh=:auto)
end  

"""
    get_bounds_for_fractions(h, fractions)

Description
Calculate the bounds containing the specified fraction(s) of area.

# Arguments:
- h:         the densities contained in the bins
- fractions: the fractional area that should be bounded

# Output:
- bounds:    the limits of the bounding area
"""
function get_bounds_for_fractions(h, fractions)
    integral = sum(h.weights)
    bounds =zeros(length(fractions))
    for (jj,fraction) in enumerate(fractions)
        minbound = 0
        maxbound = maximum(h.weights)
        newbound = 0
        for ii in 1:15
            newbound = 0.5*(minbound+maxbound)
            integral2 = sum(h.weights[h.weights.>newbound])
            newfraction = integral2/integral
            if newfraction>fraction
                minbound = newbound
            else
                maxbound = newbound
            end
        end
        bounds[jj] = newbound
    end 
    return sort(bounds)
end

"""
    create_1D_density(axis, values, range, chain_weights, fraction_1D, nbins; color, linewidth)

Description
Make the 1D density plots given the parameter values, ranges, and weights.

# Arguments:
- axis:           the axis to make the plot
- values:         the values for the x-coordinate        
- range:          the ranges for the x-coordinate        
- chain_weights:  the sample weighting from the MCMC
- fraction_1D:    the fractional area from which to compute the confidence intervals
- nbins:          number of bins, identical for all parameters   
- color:          the color of the density curve
- linewidth:      the linewidth of the density curve

# Output:
- x:              the x-coordinates of the density plot
- h:              the heights of the density plot
- y:              the normalized heights of the density plot
"""
function create_1D_density(axis, values, range, chain_weights, fraction_1D, nbins; color, linewidth)

    filter = values .> range[1] .&& values .< range[2]
    values = values[filter]
    chain_weights = weights(chain_weights[filter]) # weights is a StatsBase function
    
    h = fit(Histogram, values, chain_weights, nbins=nbins)
    x =(h.edges[1][2:end] .+ h.edges[1][1:end-1])./2
    dx = 1
    if length(x) > 1
        if x[2]-x[1] > 0
            dx = x[2]-x[1]
        end
    end
    y = h.weights/sum(h.weights*dx)
    lines!(axis, x, y, color=color, linewidth=linewidth)
    return x, h, y, dx
end

"""
    create_compound_1D_densities(axis, values_matrix, range, chain_weights_matrix, fraction_1D, nbins)

#TODO Description

# Arguments:
- axis:                 the axis to make the plot
- values_matrix:        the values for each chain of the parameter 
- range:                the ranges for the x-coordinate        
- chain_weights_matrix: the sample weighting for each chain 
- fraction_1D:          the fractional area from which to compute the confidence intervals
- nbins:                number of bins, identical for all parameters   

# Output:
- xmin:                 the left boundary of the fraction_1D area interval
- xmode:                the mode of the data, within the provided range
- xmax:                 the right boundary of the fraction_1D area interval
"""
function create_compound_1D_densities(axis, values_matrix, range, chain_weights_matrix, fraction_1D, nbins)

    # Iterate over the different chains
    for ii in 1:size(values_matrix)[1] # nchains
        values = values_matrix[ii,:]
        chain_weights = chain_weights_matrix[ii,:]
        create_1D_density(axis, values, range, chain_weights, fraction_1D, nbins, color=(:gray, 0.25), linewidth=1)
    end

    # Plot once for all the values
    x, h, y, dx = create_1D_density(axis, vec(values_matrix), range, vec(chain_weights_matrix), fraction_1D, nbins, color=(:blue, 1.0), linewidth=1)

    bound = get_bounds_for_fractions(h, [fraction_1D])[1]
    xmin = minimum(x[h.weights .>= bound]) - dx/2 # get left most value of bin
    xmax = maximum(x[h.weights .>= bound]) + dx/2
    xmode = x[argmax(h.weights)]
    if xmode - dx/2 < range[1]
        xmode = range[1]
    elseif xmode + dx/2 > range[2]
        xmode = range[2]
    end

    filter = x .>= xmin .&& x .<= xmax
    band!(axis, x[filter], zeros(length(x[filter])), y[filter], color=(:gray, 0.4))
    vlines!(axis, xmode, color=(:black, 1.0), linewidth=1)
    xlims!(axis, range[1], range[2])
    ylims!(axis, -.1*maximum(y), 5/3*maximum(y))

    return (xmin, xmode, xmax)
   
end
 
