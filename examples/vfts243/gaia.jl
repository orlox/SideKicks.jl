using Turing
using Distributions
using CairoMakie
using SideKicks


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
    d ~ Normal(50.0,30.0) # in kpc
    par_true = 1/d # in mas

    vra ~ Uniform(-500.0,500.0)
    vdec ~ Uniform(-500.0,500.0)
    pmra = vra / (4.7*d)
    pmdec = vdec / (4.7*d)
    #pmra ~ Uniform(-5.0,5.0)
    #pmdec ~ Uniform(-5.0,5.0)
    μ ~ MvNormal([par_true, pmra, pmdec],Σ)

    return (d, vra, vdec, par_true, pmra, pmdec)
end


##
nsamples = 100_000
nchains = 8
# slow step
chains = sample(sample_vfts_243_dist(μ,Σ), NUTS(0.8), MCMCThreads(), nsamples, nchains)

##
output_vals = generated_quantities(sample_vfts_243_dist(μ, Σ), chains)

##
out_vals = reduce(hcat,output_vals) # stack chains into matrix
##

all_d     = zeros(size(out_vals))
all_vra   = zeros(size(out_vals))
all_vdec  = zeros(size(out_vals))
all_par   = zeros(size(out_vals))
all_pmra  = zeros(size(out_vals))
all_pmdec = zeros(size(out_vals))

for ii in eachindex(out_vals)
    all_d[ii]     = out_vals[ii][1]
    all_vra[ii]   = out_vals[ii][2]
    all_vdec[ii]  = out_vals[ii][3]
    all_par[ii]   = out_vals[ii][4]
    all_pmra[ii]  = out_vals[ii][5]
    all_pmdec[ii] = out_vals[ii][6]
end

##

vals1 = all_par
xlabel = L"π (mas)"
vals2 = all_pmdec
ylabel = L"pm_{\delta}"
#vals2 = all_pmra
#ylabel = L"pm_{\alpha}"

chain_weights = ones(Float64, size(vals1))
fontsize=20

f = Figure()
ax = Axis(f[1,1], xlabel=xlabel, ylabel=ylabel,
          xlabelsize=fontsize, ylabelsize=fontsize)
SideKicks.create_2D_density(ax, 
                            vals1, [minimum(vals1), maximum(vals1)], 
                            vals2, [minimum(vals2), maximum(vals2)], 
                            chain_weights, 
                            [0.68,0.95,0.997],
                            120)

f
