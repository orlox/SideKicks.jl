"""
    mutable struct Observations

Observations contains the symbols, values, errors, and units of each observed parameter.
"""
@kwdef mutable struct Observations
    props::Vector{Symbol}
    vals::Vector{Float64}
    errs::Vector{Float64}
    units::Vector{Float64}
end

function Observations(obs_matrix::Vector{Vector{Any}})
    obs_matrix = stack(obs_matrix) # make into actual matrix
    return Observations(
        props = obs_matrix[1,:],
        vals  = obs_matrix[2,:],
        errs  = obs_matrix[3,:],
        units = obs_matrix[4,:])
end

macro Observations(observations)
    ex = Expr(:call)
    ex.args = [SideKicks.Observations, observations]
    return :($ex, $(string(ex)))
end