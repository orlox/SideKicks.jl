using HDF5

# TODO: close the files, but for the ExtractResults need to be able to keep the results around

function SaveResults(fname, mcmcStruct, obs_string, priors_string)
    h5open(fname, "w") do fid
        # Define results group
        results = create_group(fid, "results")
        dict_keys = keys(mcmcStruct.results)
        results["results_keys"] = string.(dict_keys)
        for key in dict_keys
            results[string(key)] = mcmcStruct.results[key]
        end
        # Define strings group
        strings  = create_group(fid, "strings")
        strings["observations"] = obs_string
        strings["priors"] = priors_string
        # Define metadata group
        meta  = create_group(fid, "metadata")
        meta["nuts_warmup_count"] = mcmcStruct.nuts_warmup_count
        meta["nuts_acceptance_rate"] = mcmcStruct.nuts_acceptance_rate
        meta["nsamples"] = mcmcStruct.nsamples
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
    return [results, observations, priors, metadata ]
end


