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
        show_CIs=true, nbins=100,
        rowcolgap=10, xticklabelrotation=pi/4,
        labelfontsize=16, tickfontsize=10, supertitlefontsize=30)

    props = plotting_props.props
    units = plotting_props.units
    names = plotting_props.names
    ranges = plotting_props.ranges

    # TODO: if extra plotting props are included that can't be used, scrap these and only plot the good ones. 
    # Print a line about this, but don't throw a warning

    # TODO: add a flag for plotting with observation or prior distribution
    # This will require making sure the props are identical

    # Confirm requested props exist
    available_props = map(Symbol, keys(results))
    for prop ∈ props
        if prop ∉ available_props
            throw(DomainError(prop, "Allowed props are only "*join([String(aprop) for aprop in available_props], ", ")))
        end
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
            axis = Axis(fig[jj+1,ii], xtickalign=1, xtickcolor = :white, ytickalign=1, ytickcolor = :white, 
                        aspect=1,
                        xlabel=names[ii], ylabel=names[jj], 
                        xlabelsize=labelfontsize, ylabelsize=labelfontsize,
                        xticklabelrotation=xticklabelrotation,
                        xticklabelsize=tickfontsize, yticklabelsize=tickfontsize)
            create_2D_density(axis, vec(results[props[ii]])/units[ii], ranges[ii], vec(results[props[jj]])/units[jj], ranges[jj], vec(results[:weights]), fractions, nbins)
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
        hideydecorations!(axis)
        if ii !=num_col
            hidexdecorations!(axis,ticks=false, minorticks=false)
        end
        str_xmode = round(xmode, sigdigits=3)
        str_upp = round(xmax-xmode, sigdigits=3)
        str_low = round(xmode-xmin, sigdigits=3)
        latex_bounds = L"%$(str_xmode)^{+%$(str_upp)}_{-%$(str_low)}"
        latex_bounds_array[ii] = latex_bounds

        println(names[ii]*"="*latex_bounds)
        if show_CIs
            # RTW: streamline this
            if ii == 1
                Label(fig[ii,ii], latex_bounds, 
                      valign=:bottom, fontsize=20,
                      tellwidth=false)

                axis = Axis(fig[ii,ii], aspect=3,
                            xgridvisible = false, ygridvisible = false)
                hidedecorations!(axis)
                hidespines!(axis)
            else
                Label(fig[ii,ii], latex_bounds, 
                      valign=:bottom, fontsize=20,
                      tellwidth=false, tellheight=false)
                axis = Axis(fig[ii,ii], aspect=1,
                            xgridvisible = false, ygridvisible = false)
                hidedecorations!(axis)
                hidespines!(axis)
            end
        end
    end     

    # Readjust rows and columns
    rowgap!(fig.layout, rowcolgap)
    colgap!(fig.layout, rowcolgap)

    if !ismissing(supertitle)
        Label(fig[0,:], text=supertitle, fontsize=supertitlefontsize)
    end

    resolution=600 .*(1,1)
    resize!(fig.scene, resolution)
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
    heatmap!(axis, y, x, h.weights, colormap=:dense)
    bounds = get_bounds_for_fractions(h, fractions)
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
    for ii in 1:size(values_matrix)[1] # nchains
        values = values_matrix[ii,:]
        chain_weights = chain_weights_matrix[ii,:]
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
 
