using Distributions

"""
    mutable struct Priors

Priors contains the prior distribution of each of the desired parameters

"""
@kwdef mutable struct Priors
    logm1_dist::Union{ContinuousUnivariateDistribution,Missing} = missing
    logm2_dist::Union{ContinuousUnivariateDistribution,Missing} = missing
    logP_dist::Union{ContinuousUnivariateDistribution,Missing} = missing 
    vkick_dist::Union{ContinuousUnivariateDistribution,Missing} = missing
    frac_dist::Union{ContinuousUnivariateDistribution,Missing} = missing 
    e_dist::Union{ContinuousUnivariateDistribution,Missing} = missing       
    venv_N_100kms_dist::Union{ContinuousUnivariateDistribution,Missing} = missing       
    venv_E_100kms_dist::Union{ContinuousUnivariateDistribution,Missing} = missing       
    venv_r_100kms_dist::Union{ContinuousUnivariateDistribution,Missing} = missing       
    rv_env_dist::Union{ContinuousUnivariateDistribution,Missing} = missing 
    pmra_dist::Union{ContinuousUnivariateDistribution,Missing} = missing  
    pmdec_dist::Union{ContinuousUnivariateDistribution,Missing} = missing 
    parallax_dist::Union{ContinuousUnivariateDistribution,Missing} = missing 
end

macro Priors(keywords...)
    ex = Expr(:call)
    for keyword in keywords
        keyword.head = :kw
    end
    ex.args = [:(SideKicks.Priors), [keyword for keyword in keywords]...]
    return :($ex, $(string(ex)))
end
