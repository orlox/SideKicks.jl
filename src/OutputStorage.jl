using HDF5

# TODO: close the files, but for the ExtractResults need to be able to keep the results around

function SaveResults(fname, kick_mcmc)
    h5open(fname, "w") do fid
        # Define results group
        results = create_group(fid, "results")
        dict_keys = keys(kick_mcmc.results)
        results["results_keys"] = string.(dict_keys)
        for key in dict_keys
            results[string(key)] = kick_mcmc.results[key]
        end
        # Define strings group
        strings  = create_group(fid, "strings")
        strings["observations"] = kick_mcmc.observations_string
        strings["priors"] = kick_mcmc.priors_string
        # Define metadata group
        meta  = create_group(fid, "metadata")
        meta["nuts_warmup_count"] = kick_mcmc.nuts_warmup_count
        meta["nuts_acceptance_rate"] = kick_mcmc.nuts_acceptance_rate
        meta["nsamples"] = kick_mcmc.nsamples
    end
end


##

function ExtractResults(fname)
    fid = h5open(fname, "r") 
    extracted_results = fid["results"]
    results = Dict()
    for key ∈ keys(extracted_results)
        results[Symbol(key)] = extracted_results[key][]
    end
    strings = fid["strings"]
    obs_string = strings["observations"][]
    observations = eval(Meta.parse(obs_string))
    priors_string = strings["priors"][]
    priors = eval(Meta.parse(priors_string))
    metadata = fid["metadata"]
    close(fid)
    return [results, observations, priors, metadata ]
end


