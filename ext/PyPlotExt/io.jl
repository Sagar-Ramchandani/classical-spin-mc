"""
--------------------------------------------------------------------------------
Common IO for all plots
--------------------------------------------------------------------------------
"""

function getFileNames(location)
	return filter((x) -> !isdir(x), readdir(location, join = true))
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
	return read.((g["mean"], g["error"]))
end

function loadObservables(fn::String, obs::Vector{Symbol})
	observables = Vector{NTuple{2, Float64}}(undef, length(obs))
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

	observables = Matrix{NTuple{2, Float64}}(undef, length(β), length(obs))
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
		groupedSpins[currentSite] = map(i -> getindex(spins, i), currentSite:nSites:totalSpins)
	end
	return groupedSpins
end