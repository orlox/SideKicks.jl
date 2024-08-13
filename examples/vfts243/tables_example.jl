using DataFrames
using CSV
using CairoMakie
using Statistics
using Turing

# this data file contains everything within 4 arcmin of VFTS 243
df = DataFrame(CSV.File(String(@__DIR__)*"/1720187991964O-result.csv"))

ruwe = df[!,"ruwe"]

pmdec = df[!,"pmdec"]
pmdec_error = df[!,"pmdec_error"]

pmra = df[!,"pmra"]
pmra_error = df[!,"pmra_error"]

parallax = df[!,"parallax"]
parallax_error = df[!,"parallax_error"]

## Apply filters
# parallax consistent with the LMC
parallax_filter = parallax.+parallax_error .> 0.02 .&& parallax.-parallax_error .< 0.02
@show sum(parallax_filter), length(parallax_filter)
# proper motions accurate to 0.1 mas/yr
accuracy_filter = pmra_error .< 0.1 .&& pmdec_error .< 0.1
@show sum(accuracy_filter), length(accuracy_filter)
# exclude runaways
non_runaway_condition = (pmra .> 1) .&& (pmra .< 2.5) .&& (pmdec .> 0) .&& (pmdec .< 2)
@show sum(non_runaway_condition), length(non_runaway_condition)
# remove bad solutions in terms of ruwe
ruwe_filter = ruwe .< 1.4

full_filter = parallax_filter .& accuracy_filter .& non_runaway_condition .& ruwe_filter
@show sum(full_filter), length(full_filter)

##

f = Figure()
ax = Axis(f[1,1])
hist!(ax, pmdec[full_filter].*(1/3600/24/365.25/1000/3600*pi/180)*(50_000*3.086e13))
@show mean(pmdec[full_filter].*(1/3600/24/365.25/1000/3600*pi/180)*(50_000*3.086e13))
f

##

f = Figure()
ax = Axis(f[1,1])
hist!(ax, pmra[full_filter].*(1/3600/24/365.25/1000/3600*pi/180)*(50_000*3.086e13))
@show mean(pmra[full_filter].*(1/3600/24/365.25/1000/3600*pi/180)*(50_000*3.086e13))
f

##
# apply the filter
pmra = pmra[full_filter]
pmra_error = pmra_error[full_filter]
pmdec = pmdec[full_filter]
pmdec_error = pmdec_error[full_filter]
par = parallax[full_filter]
par_error = parallax_error[full_filter]

par_pmra_corr = df[full_filter,"parallax_pmra_corr"]
par_pmdec_corr = df[full_filter,"parallax_pmdec_corr"]
pmra_pmdec_corr = df[full_filter,"pmra_pmdec_corr"]

##
# This model is to perform an MCMC for an observed Gaia object with
# a mean of μ and a covariance matrix Σ, using a distance prior for the LMC.
@model function sample_vfts_243_dist(μ, Σ)
    d ~ Normal(49.59,0.63) # in kpc, from https://www.nature.com/articles/s41586-019-0999-4 informally adding up systematic and statistical uncertainty
    par_true = 1/d # in mas
    pmra ~ Uniform(-5.0,5.0)
    pmdec ~ Uniform(-5.0,5.0)
    μ ~ MvNormal([par_true, pmra, pmdec],Σ)
end

# for each element in the Gaia data we will run the MCMC and extract
# the new mean and error for the velocities in the RA and DEC directions
v_ra_fix = zeros(length(pmra))
v_ra_error_fix = zeros(length(pmra))
v_dec_fix = zeros(length(pmra))
v_dec_error_fix = zeros(length(pmra))

for i in eachindex(pmra)
    μ = [parallax[i], pmra[i], pmdec[i]]
    Σ = [par_error[i]^2                                 par_pmra_corr[i]*par_error[i]*pmra_error[i]      par_pmdec_corr[i]*par_error[i]*pmdec_error[i];
        par_pmra_corr[i]*par_error[i]*pmra_error[i]     pmra_error[i]^2                                  pmra_pmdec_corr[i]*pmra_error[i]*pmdec_error[i];
        par_pmdec_corr[i]*par_error[i]*pmdec_error[i]   pmra_pmdec_corr[i]*pmra_error[i]*pmdec_error[i]  pmdec_error[i]^2]
    nsamples = 5_000
    nchains = 8
    chains = sample(sample_vfts_243_dist(μ,Σ), NUTS(0.8), MCMCThreads(), nsamples, nchains)
    # convert result to velocities in km/s
    v_ra_chain = chains[:pmra].*chains[:d].*(1/3600/24/365.25/1000/3600*pi/180).*(3.086e16)
    v_dec_chain = chains[:pmdec].*chains[:d].*(1/3600/24/365.25/1000/3600*pi/180).*(3.086e16)
    # extract new mean and sigma
    v_ra_fix[i] = mean(v_ra_chain)
    v_ra_error_fix[i] = std(v_ra_chain)
    v_dec_fix[i] = mean(v_dec_chain)
    v_dec_error_fix[i] = std(v_dec_chain)
    println("Star $i/$(length(pmra))")
    println("v_ra: $(v_ra_fix[i])+-$(v_ra_error_fix[i])")
    println("v_dec: $(v_dec_fix[i])+-$(v_dec_error_fix[i])")
end

## 

# TODO: Define here the v_r_fix and v_r_error_fix from Almeida sample

df = DataFrame(CSV.File(String(@__DIR__)*"/Almeida_truncated.csv")) 
id = df[!,"id"]
v_r_fix = df[!,"vr"]
v_r_error_fix = df[!,"vr_err"]

##
f = Figure()
ax = Axis(f[1,1])
hist!(ax, pmdec.*(1/3600/24/365.25/1000/3600*pi/180)*(50_000*3.086e13), normalization=:pdf)
hist!(ax, v_dec_fix, normalization=:pdf)
f

##
f = Figure()
ax = Axis(f[1,1])
hist!(ax, pmra.*(1/3600/24/365.25/1000/3600*pi/180)*(50_000*3.086e13), normalization = :pdf)
hist!(ax, v_ra_fix, normalization = :pdf)
f

##
# Finally, infer the distributions for the RA and DEC velocity
@model function v_model(measured_velocities, errors)
    σ_v ~ Uniform(0.0,1000.0)
    μ_v ~ Uniform(0.0,1000.0)

    for i in eachindex(measured_velocities)
        new_σ = sqrt(σ_v^2 + errors[i]^2)
        measured_velocities[i] ~ Normal(μ_v, new_σ)
    end
end

##
# first we do RA
iterations = 10_000
model = v_model(v_ra_fix, v_ra_error_fix)
chain = sample(model, NUTS(10_000,0.8), MCMCThreads(), iterations, 8);
@show median(chain[:σ_v])
@show median(chain[:μ_v])

##
# compare to RA results
f = Figure()
ax = Axis(f[1,1])
hist!(ax, v_ra_fix, normalization=:pdf)
dist = Normal(median(chain[:μ_v]), median(chain[:σ_v]))
xvals = LinRange(minimum(v_ra_fix), maximum(v_ra_fix), 100)
lines!(ax, xvals, pdf.(dist, xvals))
f

##
# and then DEC
iterations = 10_000
model = v_model(v_dec_fix, v_dec_error_fix)
chain = sample(model, NUTS(10_000,0.8), MCMCThreads(), iterations, 8);
@show median(chain[:σ_v])
@show median(chain[:μ_v])

##
# compare to DEC results
f = Figure()
ax = Axis(f[1,1])
hist!(ax, v_dec_fix, normalization=:pdf)
dist = Normal(median(chain[:μ_v]), median(chain[:σ_v]))
xvals = LinRange(minimum(v_dec_fix), maximum(v_dec_fix), 100)
lines!(ax, xvals, pdf.(dist, xvals))
f

##
# now do RV
iterations = 10_000
model = v_model(v_r_fix, v_r_error_fix)
chain = sample(model, NUTS(10_000,0.8), MCMCThreads(), iterations, 8);
@show median(chain[:σ_v])
@show median(chain[:μ_v])

##
# compare to RV results
f = Figure()
ax = Axis(f[1,1])
hist!(ax, v_r_fix, normalization=:pdf)
dist = Normal(median(chain[:μ_v]), median(chain[:σ_v]))
xvals = LinRange(minimum(v_r_fix), maximum(v_r_fix), 100)
lines!(ax, xvals, pdf.(dist, xvals))
f

##
# extract properties for VFTS 243
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
Σ = [par_err^2                          par_pmra_corr*par_err*pmra_err      par_pmdec_corr*par_err*pmdec_err;
    par_pmra_corr*par_err*pmra_err     pmra_err^2                          pmra_pmdec_corr*pmra_err*pmdec_err;
    par_pmdec_corr*par_err*pmdec_err   pmra_pmdec_corr*pmra_err*pmdec_err  pmdec_err^2]
nsamples = 5_000
nchains = 8
##
chains = sample(sample_vfts_243_dist(μ,Σ), NUTS(0.8), MCMCThreads(), nsamples, nchains)

##
sys_v_ra_chain = chains[:pmra].*chains[:d].*(1/3600/24/365.25/1000/3600*pi/180).*(3.086e16)
sys_v_dec_chain = chains[:pmdec].*chains[:d].*(1/3600/24/365.25/1000/3600*pi/180).*(3.086e16)

@show median(sys_v_ra_chain)
@show mean(sys_v_ra_chain  )
@show std(sys_v_ra_chain   )
@show median(sys_v_dec_chain)
@show mean(sys_v_dec_chain )
@show std(sys_v_dec_chain  )