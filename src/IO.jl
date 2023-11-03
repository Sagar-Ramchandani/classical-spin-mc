using HDF5
using Serialization

function writeCheckpoint!(filename::String, mc::MonteCarlo{Lattice{D,N}}) where {D,N}
    h5open(filename, "w") do f
        data = IOBuffer()
        serialize(data, mc)
        f["checkpoint"] = take!(data)
    end
end

function readCheckpoint(filename::String)
    h5open(filename, "r") do f
        data = IOBuffer(read(f["checkpoint"]))
        return deserialize(data)
    end
end

function writeUnitcell!(fn::Union{HDF5.File,HDF5.Group}, uc::UnitCell{D}) where {D}
    u = create_group(fn, "unitcell")
    #Writing primitive vectors converted to Matrix{Float64}
    u["D"] = first(typeof(uc).parameters)
    u["primitive"] = collect(reduce(hcat, uc.primitive))
    #Writing basis positions self-interactions and fields
    u["basis"] = reduce(hcat, uc.basis)
    u["interactionsField"] = reduce(hcat, uc.interactionsField)
    for i in 1:length(uc.basis)
        u["interactionsOnsite/"*string(i)] = collect(uc.interactionsOnsite[i])
    end
    #Writing two-site interactions
    #Possibly change to a single large write of flattened data
    interactions = create_group(u, "interactions")
    for i in 1:length(uc.interactions)
        interactions["$(i)/b1"] = uc.interactions[i][1][1]
        interactions["$(i)/b2"] = uc.interactions[i][1][2]
        interactions["$(i)/offset"] = collect(uc.interactions[i][2])
        interactions["$(i)/M"] = collect(uc.interactions[i][3])
    end
end

function readUnitcell(fn::Union{HDF5.File,HDF5.Group})
    u = fn["unitcell"]
    #Read primitive vectors and create unitcell
    D = read(u["D"])
    primitive = NTuple{D,SVector{D,Float64}}(eachcol(read(u["primitive"])))
    basis = Vector{SVector{D,Float64}}(eachcol(read(u["basis"])))
    interactionsField = Vector{SVector{3,Float64}}(eachcol(read(u["interactionsField"])))

    onsite = u["interactionsOnsite"]
    interactionsOnsite = [SMatrix{3,3}(read(onsite[key])) for key in keys(onsite)]

    inter = u["interactions"]

    interactions = Tuple{Pair{Int,Int},Tuple{Int,Int},SMatrix{3,3,Float64,9}
    }[(inter["$(key)/b1"] => inter["$(key)/b2"],
        inter["$(key)/offset"], inter["$(key)/M"]) for key in keys(inter)]

    return UnitCell(primitive, basis, interactions, interactionsOnsite, interactionsField)
end

function writeLattice!(fn::Union{HDF5.File,HDF5.Group}, lattice::Lattice{D,N}) where {D,N}
    """
    Only store information that cannot be reconstructed 
    with the lattice constructor
    """
    l = create_group(fn, "lattice")
    l["L"] = collect(lattice.size)
    writeUnitcell!(l, lattice.unitcell)
    l["spins"] = reduce(hcat, lattice.spins)
    return nothing
end

function readLattice(fn::Union{HDF5.File,HDF5.Group})
    """
    Reconstruct information using the lattice constructor
    """
    l = fn["lattice"]
    spins = Vector{SVector{3,Float64}}(eachcol(read(l["spins"])))
    lattice = Lattice(readUnitcell(l), Tuple(read(l["L"])))
    lattice.spins = spins
    return lattice
end

function writeMonteCarloParameters!(fn::Union{HDF5.File,HDF5.Group}, mcp::MonteCarloParameters{U}) where {U}
    p = create_group(fn, "parameters")
    #Simulation parameters
    p["beta"] = mcp.beta
    p["thermalizationSweeps"] = mcp.thermalizationSweeps
    p["measurementSweeps"] = mcp.measurementSweeps
    p["measurementRate"] = mcp.measurementRate
    p["microcanonicalRoundsPerSweep"] = mcp.microcanonicalRoundsPerSweep
    p["replicaExchangeRate"] = mcp.replicaExchangeRate
    p["randomizeInitialConfiguration"] = mcp.randomizeInitialConfiguration
    p["reportInterval"] = mcp.reportInterval
    p["checkpointInterval"] = mcp.checkpointInterval

    #Technical parameters
    p["seed"] = mcp.seed
    p["sweep"] = mcp.sweep

    p["updateFunction"] = String(Symbol(mcp.updateFunction))

    return nothing
end

function readMonteCarloParameters(fn::Union{HDF5.File,HDF5.Group})
    p = fn["parameters"]
    mcp = MonteCarloParameters(
        read(p["beta"]),
        read(p["thermalizationSweeps"]),
        read(p["measurementSweeps"]),
        read(p["measurementRate"]),
        read(p["microcanonicalRoundsPerSweep"]),
        read(p["replicaExchangeRate"]),
        read(p["randomizeInitialConfiguration"]),
        read(p["reportInterval"]),
        read(p["checkpointInterval"]),
        copy(Random.GLOBAL_RNG),
        read(p["seed"]),
        read(p["sweep"]),
        getfield(Main, Symbol(read(p["updateFunction"])))
    )
    Random.seed!(mcp.rng, mcp.seed)
    return mcp
end

function save(fn::Union{HDF5.File,HDF5.Group}, path::String, mean::Float64, error::Float64)
    fn["$(path)/mean"] = mean
    fn["$(path)/error"] = error
end

function save(fn::Union{HDF5.File,HDF5.Group}, path::String, observable::ErrorPropagator)
    save(fn, path, mean(observable[1]), std_error(observable[1]))
end

function save(fn::Union{HDF5.File,HDF5.Group}, path::String, observable::LogBinner)
    save(fn, path, mean(observable), std_error(observable))
end

function save(fn::Union{HDF5.File,HDF5.Group}, path::String, observable::FullBinner)
    save(fn, path, mean(observable), std_error(observable))
    fn["$(path)/values"] = observable.x
end

function save(fn::Union{HDF5.File,HDF5.Group}, path::String, observable::Vector{Real})
    save(fn, path, mean(observable), std_error(observable))
end

function writeObservables!(fn::Union{HDF5.File,HDF5.Group}, obs::Observables, beta::Float64, N::Float64)
    o = create_group(fn, "observables")

    for field in fieldnames(Observables)
        save(o, field, getfield(a, field))
    end

    #Save specific heat seperately
    c(e) = beta * beta * (e[2] - e[1] * e[1]) * N
    ∇c(e) = [-2.0 * beta * beta * e[1] * N, beta * beta * N]
    heat = mean(obs.energy, c)
    dheat = sqrt(abs(var(
        obs.energy, ∇c, BinningAnalysis._reliable_level(obs.energy))) /
                 obs.energy.count[BinningAnalysis._reliable_level(obs.energy)])
    save(o, "specificHeat", heat, dheat)
end

function readObservables(fn::Union{HDF5.File,HDF5.Group})
    o = fn["observables"]
    for field in fieldnames(Observables)
    end
end

function writeMonteCarlo!(filename::String, mc::MonteCarlo{Lattice{D,N}}) where {D,N}
    h5open(filename, "w") do f
        g = create_group(f, "mc")
        writeLattice!(g, mc.lattice)
        writeMonteCarloParameters!(g, mc.parameters)
        writeObservables!(g, mc.observables, mc.beta, mc.lattice.length)
    end
end

function readMonteCarlo(filename::String)
    h5open(filename, "r") do f
        g = f["mc"]
        lattice = readLattice(g)
        parameters = readMonteCarloParameters(g)
        observables = readMonteCarloObservables(g)
        return MonteCarlo(lattice, parameters, observables)
    end
end