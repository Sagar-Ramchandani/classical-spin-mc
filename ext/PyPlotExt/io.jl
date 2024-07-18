"""
--------------------------------------------------------------------------------
Common IO for all plots
--------------------------------------------------------------------------------
"""

function getFileNames(location; fileExtension = "h5")
    return filter(
        (x) -> contains(x, fileExtension) && !isdir(x), readdir(location, join = true))
end

function getFolderNames(location)
    return filter(isdir, readdir(location, join = true))
end

"""
--------------------------------------------------------------------------------
Functions for loading observables
--------------------------------------------------------------------------------
"""

function loadObservable(f::H5, obs::Symbol)
    g = f[String(obs)]
    return measurement(read.((g["mean"], g["error"]))...)
end

function loadObservables(fn::String, obs::Vector{Symbol})
    observables = Vector{Measurement{Float64}}(undef, length(obs))
    h5open(fn, "r") do f
        g = f["mc/observables"]
        for (i, o) in enumerate(obs)
            observables[i] = loadObservable(g, o)
        end
    end
    return observables
end

function loadObservables(fn::Vector{String}, obs::Vector{Symbol})
    β = map(i -> h5read(i, "mc/parameters/beta"), fn)
    perm = sortperm(β)
    β = β[perm]
    temperatures = map(inv, β)
    fn = fn[perm]

    observables = Matrix{Measurement{Float64}}(undef, length(β), length(obs))
    for (j, file) in enumerate(fn)
        h5open(file, "r") do f
            g = f["mc/observables"]
            for (i, o) in enumerate(obs)
                observables[j, i] = loadObservable(g, o)
            end
        end
    end
    return temperatures, β, observables
end

function processObservables!(temperatures::Vector{Float64}, β::Vector{Float64},
        observables, obs, f; lowT = 0.0, highT = 0.0)
    f["temperatures"] = lowT === highT ? temperatures :
                        cropXY(temperatures, temperatures, lowT, highT)[1]
    f["properties"] = String.(obs)
    g = create_group(f, "observables")
    for (i, o) in enumerate(obs)
        h = create_group(g, String(o))
        o = lowT === highT ? observables[:, i] :
            cropXY(temperatures, observables[:, i], lowT, highT)[2]
        h["val"] = getfield.(o, :val)
        h["err"] = getfield.(o, :err)
    end
end

function processObservables!(fn::Vector{String}, obs::Vector{Symbol}, saveLocation::String)
    temperatures, β, observables = loadObservables(fn, obs)
    h5open(saveLocation, "w") do f
        processObservables!(temperatures, β, observables, obs, f)
    end
end

function processObservables!(fn::Vector{Vector{String}}, obs::Vector{Symbol},
        saveLocation::String)
    data = [loadObservables(i, obs) for i in fn]
    temperatures = only(unique(getindex.(data, 1)))
    β = map(x -> 1 / x, temperatures)
    observables = getindex.(data, 3)
    meanObservables = Matrix{Measurement{Float64}}(undef, length(temperatures), length(obs))
    for i in 1:size(meanObservables, 1)
        for j in 1:size(meanObservables, 2)
            obsT = [o[i, j] for o in observables]
            meanObservables[i, j] = mean(obsT)
        end
    end
    h5open(saveLocation, "w") do f
        processObservables!(temperatures, β, meanObservables, obs, f)
    end
end

function loadProcessedObservables(fn::String)
    h5open(fn) do f
        temperatures = read(f["temperatures"])
        properties = Symbol.(read(f["properties"]))
        β = map(x -> 1 / x, temperatures)
        g = f["observables"]
        observables = Matrix{Measurement{Float64}}(
            undef, length(temperatures), length(properties))
        for (i, o) in enumerate(properties)
            h = g[String(o)]
            observables[:, i] = measurement.(read(h["val"]), read(h["err"]))
        end
        return temperatures, β, observables, properties
    end
end

"""
--------------------------------------------------------------------------------
Functions for loading spins
--------------------------------------------------------------------------------
"""

function getSpins(fn::String)
    h5open(fn, "r") do f
        s = read(f["mc/lattice/spins"])
        spins = Vector{SVector{3, Float64}}(eachcol(s))
        return spins
    end
end

function getNSites(fn::String)
    h5open(fn, "r") do f
        return size(read(f["mc/lattice/unitcell/basis"]), 2)
    end
end

function groupSpins(nSites, spins)
    groupedSpins = Vector{Vector{SVector{3, Float64}}}(undef, nSites)
    totalSpins = length(spins)
    for currentSite in 1:nSites
        groupedSpins[currentSite] = map(
            i -> getindex(spins, i), currentSite:nSites:totalSpins)
    end
    return groupedSpins
end
