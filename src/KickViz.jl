using StatsBase
using CairoMakie

export create_corner_plot, create_2D_density, create_1D_density

"""
    create_corner_plot(chain_values, names, label_names, fractions, fraction_1D, figure; 
    show_CIs=false, ranges=missing, nbins=100)

#TODO Description

# Arguments:
#TODO
- chain_values:
- names:
- label_names:
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
function create_corner_plot(chain_values, names, label_names, fractions, fraction_1D, figure; show_CIs=false, ranges=missing, nbins=100)
    ga = figure[1, 1] = GridLayout()

    if ismissing(ranges)
        ranges = [missing for i in 1:length(names)]
    end
   
    num_col = length(names)-1
    for i in 1:num_col
        for  j in i+1:num_col+1
            axis = Axis(ga[j,i],xtickalign=1,xtickcolor = :white,ytickalign=1,ytickcolor = :white, 
                   xlabel=label_names[i], ylabel=label_names[j] )
            
            create_2D_density(chain_values[names[i]], ranges[i], chain_values[names[j]], ranges[j], chain_values[:weight],fractions, axis, nbins)
            if i>1
                hideydecorations!(axis, ticks=false, minorticks=false)
            end
            if j!=num_col+1
                hidexdecorations!(axis,ticks=false, minorticks=false)
            end         
        end  
   
    end 
    for i in 1:num_col+1
        axis = Axis(ga[i,i], xgridvisible = false, ygridvisible = false,xtickalign=1,
               xlabel=label_names[i])
        (xmin, xmode, xmax) = create_1D_density(chain_values[names[i]], ranges[i], chain_values[:weight],fraction_1D,axis, nbins)
        hideydecorations!(axis)
        if i !=num_col+1
            hidexdecorations!(axis,ticks=false, minorticks=false)
        end
        if show_CIs
            axis.title = "$(xmode)^$(xmax-xmode)_$(xmode-xmin)"
        end
        print(label_names[i]*"="*"$(xmode)^$(xmax-xmode)_$(xmode-xmin)\n")
    end     
    rowgap!(ga,10)
    colgap!(ga,10)

    return figure 
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
- chain_weights:
- fractions:
- axis:
- nbins:

# Output:
#TODO
- 
"""
function create_2D_density(values1, ranges1, values2, ranges2, chain_weights, fractions, axis, nbins)
    if !ismissing(ranges1)
        filter = values1 .> ranges1[1] .&& values1 .< ranges1[2] .&&
                    values2 .> ranges2[1] .&& values2 .< ranges2[2]
        values1 = values1[filter]
        values2 = values2[filter]
        chain_weights = chain_weights[filter]
    end

    h = fit(Histogram,(values1,values2),weights(chain_weights),nbins=nbins)
    x= (h.edges[2][2:end].+h.edges[2][1:end-1])./2
    y =(h.edges[1][2:end].+h.edges[1][1:end-1])./2
    heatmap!(axis,y,x, h.weights)
    bounds = get_bounds_for_fractions(h,fractions)
   
    contour!(axis,y,x, h.weights, levels=bounds,color=:black,linewidth=2)
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
function create_1D_density(values, range, chain_weights, fraction_1D, axis, nbins)
    if !ismissing(range)
        filter = values .> range[1] .&& values .< range[2]
        values = values[filter]
        chain_weights = chain_weights[filter]
    end

    h = fit(Histogram,(values),weights(chain_weights),nbins=nbins)
    x =(h.edges[1][2:end].+h.edges[1][1:end-1])./2
    bound = get_bounds_for_fractions(h,[fraction_1D])[1]

    xmin = minimum(x[h.weights .>= bound])
    xmode = x[argmax(h.weights)]
    xmax = maximum(x[h.weights .>= bound])

    filter = x .>= xmin .&& x.<= xmax

    band!(axis, x[filter], zeros(length(x[filter])), h.weights[filter], color=(:gray,0.4))
    scatter!(axis,[xmin,xmax],[0,0])
    lines!(axis,x, h.weights)
    xlims!(axis,minimum(x),maximum(x)) 

    return (xmin, xmode, xmax)
   
end
 
