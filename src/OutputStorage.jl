using HDF5

# TODO: close the files, but for the ExtractResults need to be able to keep the results around

function SaveResults(fname, kick_mcmc)
    h5open(fname, "w") do fid
        # Define results group
        results = create_group(fid, "results")
        dict_keys = keys(kick_mcmc.results)
        results["results_keys"] = String.(dict_keys)
        for key in dict_keys
            results[String(key)] = kick_mcmc.results[key]
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

        #stats = create_group(fid, "stats")
        #stats["ess"]
        
        ess = create_group(fid, "ess")
        dict_keys = keys(kick_mcmc.ess)
        ess["ess_keys"] = String.(dict_keys)
        for key in dict_keys
            ess[String(key)] = kick_mcmc.ess[key]
        end
        
        rhat = create_group(fid, "rhat")
        dict_keys = keys(kick_mcmc.rhat)
        rhat["rhat_keys"] = String.(dict_keys)
        for key in dict_keys
            rhat[String(key)] = kick_mcmc.rhat[key]
        end
    end
end


##

function ExtractResults(fname;transpose_results=false)
    if ~isfile(fname)
        throw(ArgumentError("File not found"))
    end
    fid = h5open(fname, "r") 
    extracted_results = fid["results"]
    results = Dict()
    for key ∈ keys(extracted_results)
        results[Symbol(key)] = extracted_results[key][]
        if(transpose_results && ndims(results[Symbol(key)])==2)
            results[Symbol(key)] = transpose(results[Symbol(key)])
        end
    end
    strings = fid["strings"]
    obs_string = strings["observations"][]
    observations = eval(Meta.parse(obs_string))
    priors_string = strings["priors"][]
    priors = eval(Meta.parse(priors_string))
    metadata = fid["metadata"]

    extracted_ess = fid["ess"]
    ess = Dict()
    for key ∈ keys(extracted_ess)
        ess[Symbol(key)] = extracted_ess[key][]
    end
    extracted_rhat = fid["ess"]
    rhat = Dict()
    for key ∈ keys(extracted_rhat)
        rhat[Symbol(key)] = extracted_ess[key][]
    end
    close(fid)
    return [results, observations, priors, metadata, ess, rhat]
end


