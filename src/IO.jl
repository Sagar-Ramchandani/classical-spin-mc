using HDF5
using Serialization

#Define alias for HDF5 File or Group 
const H5 = Union{HDF5.File,HDF5.Group}

"""
--------------------------------------------------------------------------------
Saving and loading checkpoints using serialization
--------------------------------------------------------------------------------
"""

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

"""
--------------------------------------------------------------------------------
Saving and loading unitcells
--------------------------------------------------------------------------------
"""

function writeUnitcell!(fn::H5, uc::UnitCell{D}) where {D}
    u = create_group(fn, "unitcell")
    #Writing primitive vectors converted to Matrix{Float64}
    u["D"] = dimension(uc)
    u["primitive"] = collect(reduce(hcat, uc.primitive))
    #Writing basis positions self-interactions and fields
    u["basis"] = reduce(hcat, uc.basis)
    u["interactionsField"] = reduce(hcat, uc.interactionsField)
    for i in 1:length(uc.basis)
        u["interactionsOnsite/"*string(i)] = collect(uc.interactionsOnsite[i])
    end
    #Writing two-site interactions
    interactions = create_group(u, "interactions")
    for i in 1:length(uc.interactions)
        interactions["$(i)/b1"] = uc.interactions[i][1][1]
        interactions["$(i)/b2"] = uc.interactions[i][1][2]
        interactions["$(i)/offset"] = collect(uc.interactions[i][2])
        interactions["$(i)/M"] = collect(uc.interactions[i][3])
    end

    u["anisotropyFunction"] = String(Symbol(uc.anisotropyFunction))
    u["anisotropyParameters"] = uc.anisotropyParameters

    return nothing
end

function readUnitcell(fn::H5)
    u = fn["unitcell"]
    #Read primitive vectors and create unitcell
    D = read(u["D"])
    primitive = NTuple{D,SVector{D,Float64}}(eachcol(read(u["primitive"])))
    basis = Vector{SVector{D,Float64}}(eachcol(read(u["basis"])))
    interactionsField = Vector{SVector{3,Float64}}(eachcol(read(u["interactionsField"])))

    onsite = u["interactionsOnsite"]
    interactionsOnsite = [SMatrix{3,3}(read(onsite[key])) for key in keys(onsite)]

    inter = u["interactions"]

    interactions = [(read(inter["$(key)/b1"]) => read(inter["$(key)/b2"]),
        NTuple{D,Int}(read(inter["$(key)/offset"])),
        SMatrix{3,3,Float64,9}(read(inter["$(key)/M"])))
                    for key in keys(inter)]

    anisotropyFunction = getfield(Main, Symbol(read(u["anisotropyFunction"])))
    anisotropyParameters = read(u["anisotropyParameters"])

    return UnitCell(primitive, basis, interactions, interactionsOnsite, interactionsField,
        anisotropyFunction, anisotropyParameters)
end

"""
--------------------------------------------------------------------------------
Saving and loading lattices
--------------------------------------------------------------------------------
"""

function writeLattice!(fn::H5, lattice::Lattice{D,N}) where {D,N}
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

function readLattice(fn::H5)
    """
    Reconstruct information using the lattice constructor
    """
    l = fn["lattice"]
    spins = Vector{SVector{3,Float64}}(eachcol(read(l["spins"])))
    lattice = Lattice(readUnitcell(l), Tuple(read(l["L"])))
    lattice.spins = spins
    return lattice
end

"""
--------------------------------------------------------------------------------
Saving and loading MonteCarloParameters
--------------------------------------------------------------------------------
"""

function writeMonteCarloParameters!(fn::H5, mcp::MonteCarloParameters{U}) where {U}
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

    #Saving state of a general rng
    rng = create_group(p, "rng")
    rng["type"] = String(Symbol(typeof(mcp.rng)))
    currentRNG = copy(mcp.rng)
    for field in fieldnames(typeof(currentRNG))
        rng[String(field)] = getfield(mcp.rng, field)
    end

    #Technical parameters
    p["seed"] = mcp.seed
    p["sweep"] = mcp.sweep

    p["updateFunction"] = String(Symbol(mcp.updateFunction))
    p["updateParameter"] = mcp.updateParameter

    return nothing
end

function readMonteCarloParameters(fn::H5)
    p = fn["parameters"]
    rng = p["rng"]

    """
    Note: The code needs the rng type and also the update function
    to be available in the Main namespace. 
    """
    #Load the state of a general rng 
    RandomGenerator = getfield(Main, Symbol(read(rng["type"])))()
    for field in fieldnames(typeof(RandomGenerator))
        setfield!(RandomGenerator, field, read(rng[String(field)]))
    end
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
        RandomGenerator,
        read(p["seed"]),
        read(p["sweep"]),
        getfield(Main, Symbol(read(p["updateFunction"]))),
        read(p["updateParameter"])
    )
    return mcp
end

"""
--------------------------------------------------------------------------------
Saving and loading MonteCarloStatistics
--------------------------------------------------------------------------------
"""

function writeMonteCarloStatistics!(fn::H5, mcs::MonteCarloStatistics)
    s = create_group(fn, "statistics")
    fields = fieldnames(typeof(mcs))
    for field in fields
        s[String(field)] = getfield(mcs, field)
    end
    return nothing
end

function readMonteCarloStatistics(fn::H5)
    s = fn["statistics"]
    fields = fieldnames(MonteCarloStatistics)
    return MonteCarloStatistics([read(s[String(field)]) for field in fields]...)
end

"""
--------------------------------------------------------------------------------
save! functions for observable fields
--------------------------------------------------------------------------------
"""

#Saving general observable with mean and error
function save!(fn::H5, path::String, mean::T, error::T) where {T}
    fn["$(path)/mean"] = mean
    fn["$(path)/error"] = error
    return nothing
end

function save!(fn::H5, path::String, meanerror::Tuple{Float64,Float64})
    save!(fn, path, meanerror...)
end

#Saving accumulator
function save!(fn::H5, path::String, accumulators::NTuple{D,BinningAnalysis.Variance{T}}) where {D,T}
    for (i, accum) in enumerate(accumulators)
        acc = "$(path)/accumulators/$(i)/"
        fn[acc*"δ"] = accum.δ
        fn[acc*"m1"] = accum.m1
        fn[acc*"m2"] = accum.m2
        fn[acc*"count"] = accum.count
    end
    return nothing
end

#Saving compressor
function save!(fn::H5, path::String, compressors::NTuple{D,BinningAnalysis.Compressor{T}}) where {D,T}
    for (i, comp) in enumerate(compressors)
        fn["$(path)/compressors/$(i)/switch"] = comp.switch
        fn["$(path)/compressors/$(i)/value"] = comp.value
    end
    return nothing
end

#Saving EPCompressor
function save!(fn::H5, path::String, compressors::NTuple{D,BinningAnalysis.EPCompressor{T}}) where {D,T}
    for (i, comp) in enumerate(compressors)
        fn["$(path)/compressors/$(i)/switch"] = comp.switch
        fn["$(path)/compressors/$(i)/values"] = comp.values
    end
    return nothing
end

#Saving ErrorPropagator
function save!(fn::H5, path::String, observable::ErrorPropagator{T,D}) where {T,D}
    save!(fn, path, means(observable)[1], std_errors(observable)[1])
    save!(fn, path, observable.compressors)
    fn["$(path)/sums1D"] = reduce(hcat, observable.sums1D)
    fn["$(path)/sums2D"] = reduce(hcat, observable.sums2D)
    fn["$(path)/count"] = observable.count
    return nothing
end

#Saving LogBinner
function save!(fn::Union{HDF5.File,HDF5.Group}, path::String, observable::LogBinner)
    save!(fn, path, mean(observable), std_error(observable))
    save!(fn, path, observable.accumulators)
    save!(fn, path, observable.compressors)
    return nothing
end

#Saving FullBinner
function save!(fn::Union{HDF5.File,HDF5.Group}, path::String, observable::FullBinner)
    save!(fn, path, mean(observable), std_error(observable))
    fn["$(path)/values"] = reduce(hcat, observable.x)
    return nothing
end

#Saving Vector
function save!(fn::Union{HDF5.File,HDF5.Group}, path::String, observable::Vector{T}) where {T<:Real}
    save!(fn, path, mean(observable), std(observable))
    fn["$(path)/values"] = observable
    return nothing
end

"""
--------------------------------------------------------------------------------
load functions for observable fields
--------------------------------------------------------------------------------
"""

#Loading a general observable with mean and error
function load(::Val{T}, fn::H5, path::String) where {T}
    mean = read(fn["$(path)/mean"])
    error = read(fn["$(path)/error"])
    return (mean, error)
end

#Loading a vector
function load(::Val{Vector{T}}, fn::H5, path::String) where {T<:Real}
    return read(fn["$(path)/values"])
end

#Loading accumulator
function load(::Val{BinningAnalysis.Variance{T}}, fn::H5, path::String) where {T}
    indexes = map((x) -> parse(Int, x), keys(fn["$(path)/accumulators"]))
    Nindexes = length(indexes)
    accums = Vector{BinningAnalysis.Variance{T}}(undef, Nindexes)
    for i in indexes
        acc = fn["$(path)/accumulators/$(i)/"]
        accums[i] = BinningAnalysis.Variance(read(acc["δ"]), read(acc["m1"]), read(acc["m2"]), read(acc["count"]))
    end
    return NTuple{Nindexes,BinningAnalysis.Variance{T}}(accums)
end

#Loading compressor
function load(::Val{BinningAnalysis.Compressor{T}}, fn::H5, path::String) where {T}
    indexes = map((x) -> parse(Int, x), keys(fn["$(path)/compressors"]))
    Nindexes = length(indexes)
    comps = Vector{BinningAnalysis.Compressor{T}}(undef, Nindexes)
    for i in indexes
        cp = fn["$(path)/compressors/$(i)/"]
        comps[i] = BinningAnalysis.Compressor(read(cp["value"]), read(cp["switch"]))
    end
    return NTuple{Nindexes,BinningAnalysis.Compressor{T}}(comps)
end

#Loading EPCompressor
function load(::Val{BinningAnalysis.EPCompressor{T}}, fn::H5, path::String) where {T}
    indexes = map((x) -> parse(Int, x), keys(fn["$(path)/compressors"]))
    Nindexes = length(indexes)
    comps = Vector{BinningAnalysis.EPCompressor{T}}(undef, Nindexes)
    for i in indexes
        cp = fn["$(path)/compressors/$(i)/"]
        comps[i] = BinningAnalysis.EPCompressor(read(cp["values"]), read(cp["switch"]))
    end
    return NTuple{Nindexes,BinningAnalysis.EPCompressor{T}}(comps)
end

#Loading ErrorPropagator
function load(::Val{ErrorPropagator{T,D}}, fn::H5, path::String) where {T,D}
    comp = load(Val(BinningAnalysis.EPCompressor{T}), fn, path)
    sums1D = eachcol(read(fn["$(path)/sums1D"]))
    sums1D = NTuple{length(sums1D),Vector{T}}(sums1D)
    sums2D = (read(fn["$(path)/sums2D"]))
    dimensions = size(sums2D, 1)
    sums2D = reshape(sums2D, dimensions, dimensions, :)
    dimensions = size(sums2D, 3)
    sums2D = NTuple{dimensions,Matrix{T}}([Matrix{T}(sums2D[:, :, i]) for i in 1:dimensions])
    count = read(fn["$(path)/count"])
    return ErrorPropagator(comp, sums1D, sums2D, count)
end

#Loading LogBinner
function load(::Val{LogBinner{T,D,V}}, fn::Union{HDF5.File,HDF5.Group}, path::String) where {T,D,V}
    accum = load(Val(BinningAnalysis.Variance{T}), fn, path)
    comp = load(Val(BinningAnalysis.Compressor{T}), fn, path)
    return LogBinner{T,D}(comp, accum)
end

#Loading FullBinner
function load(::Val{FullBinner{T}}, fn::Union{HDF5.File,HDF5.Group}, path::String) where {T}
    return FullBinner(T(ncols(read(fn["$(path)/values"]))))
end

#Top-level load function
function load(fn::H5, path::String)
    observableType = eval(Meta.parse(read(fn["$(path)/observableType"])))
    return load(Val(observableType), fn, path)
end

"""
Functions below are experimental
"""

"""
--------------------------------------------------------------------------------
Saving and loading Observables
--------------------------------------------------------------------------------
"""

function writeObservables!(fn::Union{HDF5.File,HDF5.Group}, obs::Observables)
    o = create_group(fn, "observables")

    for field in fieldnames(Observables)
        currentObservableName = String(field)
        currentObservable = getfield(obs, field)
        save!(o, String(field), getfield(obs, field))
        o["$(currentObservableName)/observableType"] = String(Symbol(typeof(currentObservable)))
    end

    return nothing
end

function readObservables(fn::Union{HDF5.File,HDF5.Group})
    o = fn["observables"]
    return Observables([load(o, String(field)) for field in fieldnames(Observables)]...)
end

"""
--------------------------------------------------------------------------------
Saving and loading MonteCarlo
--------------------------------------------------------------------------------
"""

function writeMonteCarlo!(fn::H5, mc::MonteCarlo)
    g = create_group(fn, "mc")
    writeLattice!(g, mc.lattice)
    writeMonteCarloParameters!(g, mc.parameters)
    writeObservables!(g, mc.observables)
    return nothing
end

function writeMonteCarlo!(filename::String, mc::MonteCarlo)
    h5open(filename, "w") do f
        writeMonteCarlo!(f, mc)
    end
    return nothing
end

function readMonteCarlo(fn::H5)
    g = fn["mc"]
    lattice = readLattice(g)
    parameters = readMonteCarloParameters(g)
    observables = readObservables(g)
    return MonteCarlo(lattice, parameters, observables)
end

function readMonteCarlo(filename::String)
    h5open(filename, "r") do f
        return readMonteCarlo(f)
    end
end