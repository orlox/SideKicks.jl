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


"""
    create_corner_plot(chain_values, names, names, fractions, fraction_1D, figure; 
    show_CIs=false, ranges=missing, nbins=100)

#TODO Description

# Arguments:
#TODO
- chain_values:
- names:
- names:
- fractions:
- fraction_1D:
- figure:
- show_CIs=false:
- ranges=missing:
- nbins=100:

# Output:
#TODO
- 
"""

function create_corner_plot(results, plotting_props; 
        fig=Figure(), supertitle=missing,
        fractions=[0.68,0.95,0.997], fraction_1D=0.68, 
        show_CIs=false, nbins=100,
        rowcolgap=15, xticklabelrotation=20,
        labelfontsize=16, tickfontsize=4, supertitlefontsize=30)

    props = plotting_props.props
    units = plotting_props.units
    names = plotting_props.names
    ranges = plotting_props.ranges

    # Confirm requested props exist
    available_props = keys(results)
    for prop ∈ props
        if prop ∉ available_props
            throw(DomainError(prop, "Allowed props are only "*join([String(aprop) for aprop in available_props], ", ")))
        end
    end

    # Add ranges if none supplied
    num_props = length(props)
    for i in 1:num_props
        values = vec(results[props[i]])/units[i]
        minval = minimum(values)
        maxval = maximum(values)
        if ismissing(ranges[i]) || minval > ranges[i][2] || maxval < ranges[i][1]
            ranges[i] = [minval, maxval]
            print("Range set for ", props[i], ": ", ranges[i])
        end
    end

    # Create 2D density plots
    num_col = num_props
    for i in 1:num_col-1
        for  j in i+1:num_col
            axis = Axis(fig[j,i], xtickalign=1, xtickcolor = :white, ytickalign=1, ytickcolor = :white, 
                        xlabel=names[i], ylabel=names[j], 
                        xlabelsize=labelfontsize, ylabelsize=labelfontsize,
                        xticklabelrotation=xticklabelrotation,
                        xticklabelsize=tickfontsize, yticklabelsize=tickfontsize)
            create_2D_density(axis, vec(results[props[i]])/units[i], ranges[i], vec(results[props[j]])/units[j], ranges[j], vec(results[:weights]), fractions, nbins)
            if i>1
                hideydecorations!(axis, ticks=false, minorticks=false)
            end
            if j!=num_col
                hidexdecorations!(axis,ticks=false, minorticks=false)
            end         
        end  
    end 
    # Create 1D PDFs along the diagonal
    for i in 1:num_col
        axis = Axis(fig[i,i], xgridvisible = false, ygridvisible = false, xtickalign=1, xlabel=names[i],
                        xlabelsize=labelfontsize, ylabelsize=labelfontsize,
                        xticklabelrotation=xticklabelrotation,
                        xticklabelsize=tickfontsize, yticklabelsize=tickfontsize)
        (xmin, xmode, xmax) = create_compound_1D_densities(axis, results[props[i]]/units[i], ranges[i], results[:weights], fraction_1D, nbins)
        hideydecorations!(axis)
        if i !=num_col
            hidexdecorations!(axis,ticks=false, minorticks=false)
        end
        str_xmode = round(xmode, sigdigits=3)
        str_upp = round(xmax-xmode, sigdigits=3)
        str_low = round(xmode-xmin, sigdigits=3)
        latex_bounds = L"%$(str_xmode)^{+%$(str_upp)}_{-%$(str_low)}"
        if show_CIs
            axis.title = latex_bounds
        end
        println(names[i]*"="*latex_bounds)
    end     
    rowgap!(fig.layout, rowcolgap)
    colgap!(fig.layout, rowcolgap)

    if !ismissing(supertitle)
        Label(fig[0,:], text=supertitle, fontsize=supertitlefontsize)
    end

    return fig
end

"""
    create_2D_density(values1, ranges1, values2, ranges2, chain_weights, fractions, axis, nbins)

#TODO Description

# Arguments:
#TODO
- values1:
- ranges1:
- values2:
- ranges2:
- weights:
- fractions:
- axis:
- nbins:

# Output:
#TODO
- 
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
    heatmap!(axis, y, x, h.weights)
    bounds = get_bounds_for_fractions(h, fractions)
   
    contour!(axis, y, x, h.weights, levels=bounds, color=:black, linewidth=2)
end  

"""
    get_bounds_for_fractions(h, fractions)

#TODO Description

# Arguments:
#TODO
h:
fractions:

# Output:
#TODO
- 
"""
function get_bounds_for_fractions(h, fractions)
    integral = sum(h.weights)
    bounds =zeros(length(fractions))
    for (j,fraction) in enumerate(fractions)
        minbound = 0
        maxbound = maximum(h.weights)
        newbound = 0
        for i in 1:15
            newbound = 0.5*(minbound+maxbound)
            integral2 = sum(h.weights[h.weights.>newbound])
            newfraction = integral2/integral
            if newfraction>fraction
                minbound = newbound
            else
                maxbound = newbound
            end

        end
        bounds[j] = newbound
        
    end 
    return bounds
end

"""
    create_1D_density(values, range, chain_weights, fraction_1D, axis, nbins)

#TODO Description

# Arguments:
#TODO
- values:
- range:
- chain_weights:
- fraction_1D:
- axis:
- nbins:

# Output:
#TODO
- 
"""
function create_1D_density(axis, values, range, chain_weights, fraction_1D, nbins; color, linewidth)

    filter = values .> range[1] .&& values .< range[2]
    values = values[filter]
    chain_weights = weights(chain_weights[filter]) # weights is a StatsBase function
    
    h = fit(Histogram, values, chain_weights, nbins=nbins)
    x =(h.edges[1][2:end] .+ h.edges[1][1:end-1])./2
    lines!(axis, x, h.weights/sum(h.weights), color=color, linewidth=linewidth)
    return x, h
end

"""
    create_1D_density(values, range, chain_weights, fraction_1D, axis, nbins)

#TODO Description

# Arguments:
#TODO
- values:
- range:
- chain_weights:
- fraction_1D:
- axis:
- nbins:

# Output:
#TODO
- 
"""
function create_compound_1D_densities(axis, values_matrix, range, chain_weights_matrix, fraction_1D, nbins)

    # Iterate over the different chains
    for i in 1:size(values_matrix)[1] # nchains
        values = values_matrix[i,:]
        chain_weights = chain_weights_matrix[i,:]
        create_1D_density(axis, values, range, chain_weights, fraction_1D, nbins, color=(:gray, 0.4), linewidth=1)
    end

    # Plot once for all the values
    x, h = create_1D_density(axis, vec(values_matrix), range, vec(chain_weights_matrix), fraction_1D, nbins, color=(:blue, 1.0), linewidth=2)

    bound = get_bounds_for_fractions(h, [fraction_1D])[1]
    xmin = minimum(x[h.weights .>= bound])
    xmode = x[argmax(h.weights)]
    xmax = maximum(x[h.weights .>= bound])

    filter = x .>= xmin .&& x.<= xmax
    band!(axis, x[filter], zeros(length(x[filter])), h.weights[filter]/sum(h.weights), color=(:gray, 0.4))
    lines!(axis, x, h.weights/sum(h.weights))
    vlines!(axis, xmode, color=(:black, 1.0), linewidth=0.5)
    xlims!(axis, range[1], range[2])

    return (xmin, xmode, xmax)
   
end
 
