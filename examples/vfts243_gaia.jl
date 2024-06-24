using Turing
using Distributions
using CairoMakie

par = -0.04681187690491371
pmra = 1.7220500229370916
pmdec = 0.6031479075324057

par_err = 0.02444555
pmra_err = 0.033650085
pmdec_err = 0.031851716

par_pmra_corr = 0.19750258
par_pmdec_corr = -0.1501674
pmra_pmdec_corr = -0.0039443113

μ = [par, pmra, pmdec]
Σ = [par_err^2                         par_pmra_corr*par_err*pmra_err      par_pmdec_corr*par_err*pmdec_err;
    par_pmra_corr*par_err*pmra_err     pmra_err^2                          pmra_pmdec_corr*pmra_err*pmdec_err;
    par_pmdec_corr*par_err*pmdec_err   pmra_pmdec_corr*pmra_err*pmdec_err  pmdec_err^2]

@model function sample_vfts_243_dist(μ, Σ)
    x ~ MvNormal(μ,Σ)
end

##
nsamples = 100_000
chains = sample(sample_vfts_243_dist(μ,Σ), NUTS(0.8), 100_000)
##

chain_weights = ones(Float64, 1, nsamples)
f = Figure()
ax = Axis(f[1,1])
SideKicks.create_2D_density(ax, chains["x[1]"], [], chains["x[2]"], [], 
                            chain_weights, nbins=120)
