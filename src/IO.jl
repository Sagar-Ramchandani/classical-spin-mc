"""
    const H5 = Union{HDF5.File,HDF5.Group}
Used to define an alias for a HDF5 File or Group.
"""
const H5 = Union{HDF5.File, HDF5.Group}

"""
--------------------------------------------------------------------------------
Saving and loading checkpoints using serialization
--------------------------------------------------------------------------------
"""

"""
    function writeCheckpoint!(filename::String, mc::MonteCarlo{Lattice{D,N}}) where {D,N}
Interface function to write a MonteCarlo checkpoint.

!!! warning "Checkpoints"
    Checkpoints use Serialization and thus are not recommended for long term storage 
    as they may not be supported on a different Julia version.
"""
function writeCheckpoint!(filename::String,
        mc::MonteCarlo{T, P, O}) where {
        T <: Lattice, P <: MonteCarloParameters, O <: AbstractObservables}
    h5open(filename, "w") do f
        data = IOBuffer()
        serialize(data, mc)
        f["checkpoint"] = take!(data)
    end
end

"""
    function readCheckpoint(filename::String)
Interface function to read a MonteCarlo checkpoint.

!!! warning "Checkpoints"
    Checkpoints use Serialization and thus are not recommended for long term storage 
    as they may not be supported on a different Julia version.
"""
function readCheckpoint(filename::String)
    h5open(filename, "r") do f
        data = IOBuffer(read_dataset(f, "checkpoint"))
        return deserialize(data)
    end
end

"""
--------------------------------------------------------------------------------
Saving and loading unitcells
--------------------------------------------------------------------------------
"""

"""
    function writeUnitcell!(fn::H5, uc::UnitCell{D}) where {D}
Writes the unitcell in a H5 object.
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
        u["interactionsOnsite/" * string(i)] = collect(uc.interactionsOnsite[i])
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

"""
    function readUnitcell(fn::H5)
Reads the unitcell from a H5 object.
"""
function readUnitcell(fn::H5)
    u = open_group(fn, "unitcell")
    #Read primitive vectors and create unitcell
    D = read_dataset(u, "D")
    primitive = NTuple{D, SVector{D, Float64}}(eachcol(read_dataset(u, "primitive")))
    basis = Vector{SVector{D, Float64}}(eachcol(read_dataset(u, "basis")))
    interactionsField = Vector{SVector{3, Float64}}(eachcol(read_dataset(
        u, "interactionsField")))

    onsite = open_group(u, "interactionsOnsite")
    interactionsOnsite = [SMatrix{3, 3}(read_dataset(onsite, key)) for key in keys(onsite)]

    inter = open_group(u, "interactions")

    interactions = [(read_dataset(inter, "$(key)/b1") => read_dataset(inter, "$(key)/b2"),
                        NTuple{D, Int}(read_dataset(inter, "$(key)/offset")),
                        SMatrix{3, 3, Float64, 9}(read_dataset(inter, "$(key)/M")))
                    for key in keys(inter)]

    anisotropyFunction = getfield(Main, Symbol(read_dataset(u, "anisotropyFunction")))
    anisotropyParameters = read_dataset(u, "anisotropyParameters")

    return UnitCell(primitive, basis, interactions, interactionsOnsite, interactionsField,
        anisotropyFunction, anisotropyParameters)
end

"""
--------------------------------------------------------------------------------
Saving and loading lattices
--------------------------------------------------------------------------------
"""

"""
    function writeLattice!(fn::H5, lattice::Lattice{D,N}) where {D,N}
Writes the Lattice in a H5 object.
"""
function writeLattice!(fn::H5, lattice::Lattice{D, N}) where {D, N}
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

"""
    function readLattice(fn::H5)
Reads the Lattice from a H5 object.
"""
function readLattice(fn::H5)
    """
    Reconstruct information using the lattice constructor
    """
    l = open_group(fn, "lattice")
    spins = Vector{SVector{3, Float64}}(eachcol(read_dataset(l, "spins")))
    lattice = Lattice(readUnitcell(l), Tuple(read_dataset(l, "L")))
    lattice.spins = spins
    return lattice
end

"""
--------------------------------------------------------------------------------
Saving and loading MonteCarloParameters
--------------------------------------------------------------------------------
"""

"""
    function writeMonteCarloParameters!(fn::H5, mcp::MonteCarloParameters{U}) where {U}
Writes the MonteCarloParameters in a H5 object.
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

"""
    function readMonteCarloParameters(fn::H5)
Reads the MonteCarloParameters from a H5 object.
"""
function readMonteCarloParameters(fn::H5)
    p = open_group(fn, "parameters")
    rng = open_group(p, "rng")

    """
    Note: The code needs the rng type and also the update function
    to be available in the Main namespace. 
    """
    #Load the state of a general rng 
    #RandomGenerator = getfield(Main, Symbol(read(rng["type"])))()
    RandomGenerator = getfield(Main, Symbol(last(split(read_dataset(rng, "type"), '.'))))()
    for field in fieldnames(typeof(RandomGenerator))
        setfield!(RandomGenerator, field, read_dataset(rng, String(field)))
    end
    mcp = MonteCarloParameters(
        read_dataset(p, "beta"),
        read_dataset(p, "thermalizationSweeps"),
        read_dataset(p, "measurementSweeps"),
        read_dataset(p, "measurementRate"),
        read_dataset(p, "microcanonicalRoundsPerSweep"),
        read_dataset(p, "replicaExchangeRate"),
        read_dataset(p, "randomizeInitialConfiguration"),
        read_dataset(p, "reportInterval"),
        read_dataset(p, "checkpointInterval"),
        RandomGenerator,
        read_dataset(p, "seed"),
        read_dataset(p, "sweep"),
        getfield(Main, Symbol(read_dataset(p, "updateFunction"))),
        read_dataset(p, "updateParameter")
    )
    return mcp
end

"""
--------------------------------------------------------------------------------
Saving and loading MonteCarloStatistics
--------------------------------------------------------------------------------
"""

"""
    function writeMonteCarloStatistics!(fn::H5, mcs::MonteCarloStatistics)
Writes the MonteCarloStatistics in a H5 object.
"""
function writeMonteCarloStatistics!(fn::H5, mcs::MonteCarloStatistics)
    s = create_group(fn, "statistics")
    fields = fieldnames(typeof(mcs))
    for field in fields
        s[String(field)] = getfield(mcs, field)
    end
    return nothing
end

"""
    function readMonteCarloStatistics(fn::H5)
Reads the MonteCarloStatistics from a H5 object.
"""
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

function save!(fn::H5, path::String, meanerror::Tuple{Float64, Float64})
    save!(fn, path, meanerror...)
end

#Saving accumulator
function save!(fn::H5, path::String,
        accumulators::NTuple{D, BinningAnalysis.Variance{T}}) where {D, T}
    for (i, accum) in enumerate(accumulators)
        acc = "$(path)/accumulators/$(i)/"
        fn[acc * "δ"] = accum.δ
        fn[acc * "m1"] = accum.m1
        fn[acc * "m2"] = accum.m2
        fn[acc * "count"] = accum.count
    end
    return nothing
end

#Saving compressor
function save!(fn::H5, path::String,
        compressors::NTuple{D, BinningAnalysis.Compressor{T}}) where {D, T}
    for (i, comp) in enumerate(compressors)
        fn["$(path)/compressors/$(i)/switch"] = comp.switch
        fn["$(path)/compressors/$(i)/value"] = comp.value
    end
    return nothing
end

#Saving EPCompressor
function save!(fn::H5, path::String,
        compressors::NTuple{D, BinningAnalysis.EPCompressor{T}}) where {D, T}
    for (i, comp) in enumerate(compressors)
        fn["$(path)/compressors/$(i)/switch"] = comp.switch
        fn["$(path)/compressors/$(i)/values"] = comp.values
    end
    return nothing
end

#Saving ErrorPropagator
function save!(fn::H5, path::String, observable::ErrorPropagator{T, D}) where {T, D}
    save!(fn, path, means(observable)[1], std_errors(observable)[1])
    save!(fn, path, observable.compressors)
    fn["$(path)/sums1D"] = reduce(hcat, observable.sums1D)
    fn["$(path)/sums2D"] = reduce(hcat, observable.sums2D)
    fn["$(path)/count"] = observable.count
    return nothing
end

#Saving LogBinner
function save!(fn::Union{HDF5.File, HDF5.Group}, path::String, observable::LogBinner)
    save!(fn, path, mean(observable), std_error(observable))
    save!(fn, path, observable.accumulators)
    save!(fn, path, observable.compressors)
    return nothing
end

#Saving FullBinner
function save!(fn::Union{HDF5.File, HDF5.Group}, path::String, observable::FullBinner)
    save!(fn, path, mean(observable), std_error(observable))
    fn["$(path)/values"] = reduce(hcat, observable.x)
    return nothing
end

#Saving Vector
function save!(fn::Union{HDF5.File, HDF5.Group}, path::String,
        observable::Vector{T}) where {T <: Real}
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
    mean = read_dataset(fn, "$(path)/mean")
    error = read_dataset(fn, "$(path)/error")
    return (mean, error)
end

#Loading a vector
function load(::Val{Vector{T}}, fn::H5, path::String) where {T <: Real}
    return read_dataset(fn, "$(path)/values")
end

#Loading accumulator
function load(::Val{BinningAnalysis.Variance{T}}, fn::H5, path::String) where {T}
    indexes = map((x) -> parse(Int, x), keys(open_group(fn, "$(path)/accumulators")))
    Nindexes = length(indexes)
    accums = Vector{BinningAnalysis.Variance{T}}(undef, Nindexes)
    for i in indexes
        acc = open_group(fn, "$(path)/accumulators/$(i)/")
        accums[i] = BinningAnalysis.Variance(
            read_dataset(acc, "δ"),
            read_dataset(acc, "m1"),
            read_dataset(acc, "m2"),
            read_dataset(acc, "count"))
    end
    return NTuple{Nindexes, BinningAnalysis.Variance{T}}(accums)
end

#Loading compressor
function load(::Val{BinningAnalysis.Compressor{T}}, fn::H5, path::String) where {T}
    indexes = map((x) -> parse(Int, x), keys(open_group(fn, "$(path)/compressors")))
    Nindexes = length(indexes)
    comps = Vector{BinningAnalysis.Compressor{T}}(undef, Nindexes)
    for i in indexes
        cp = open_group(fn, "$(path)/compressors/$(i)/")
        comps[i] = BinningAnalysis.Compressor(
            read_dataset(cp, "value"), read_dataset(cp, "switch"))
    end
    return NTuple{Nindexes, BinningAnalysis.Compressor{T}}(comps)
end

#Loading EPCompressor
function load(::Val{BinningAnalysis.EPCompressor{T}}, fn::H5, path::String) where {T}
    indexes = map((x) -> parse(Int, x), keys(open_group(fn, "$(path)/compressors")))
    Nindexes = length(indexes)
    comps = Vector{BinningAnalysis.EPCompressor{T}}(undef, Nindexes)
    for i in indexes
        cp = open_group(fn, "$(path)/compressors/$(i)/")
        comps[i] = BinningAnalysis.EPCompressor(
            read_dataset(cp, "values"), read_dataset(cp, "switch"))
    end
    return NTuple{Nindexes, BinningAnalysis.EPCompressor{T}}(comps)
end

#Loading ErrorPropagator
function load(::Val{ErrorPropagator{T, D}}, fn::H5, path::String) where {T, D}
    comp = load(Val(BinningAnalysis.EPCompressor{T}), fn, path)
    sums1D = eachcol(read_dataset(fn, "$(path)/sums1D"))
    sums1D = NTuple{length(sums1D), Vector{T}}(sums1D)
    sums2D = read_dataset(fn, "$(path)/sums2D")
    dimensions = size(sums2D, 1)
    sums2D = reshape(sums2D, dimensions, dimensions, :)
    dimensions = size(sums2D, 3)
    sums2D = NTuple{dimensions, Matrix{T}}([Matrix{T}(sums2D[:, :, i])
                                            for i in 1:dimensions])
    count = read_dataset(fn, "$(path)/count")
    return ErrorPropagator(comp, sums1D, sums2D, count)
end

#Loading LogBinner
function load(::Val{LogBinner{T, D, V}}, fn::Union{HDF5.File, HDF5.Group},
        path::String) where {T, D, V}
    accum = load(Val(BinningAnalysis.Variance{T}), fn, path)
    comp = load(Val(BinningAnalysis.Compressor{T}), fn, path)
    return LogBinner{T, D}(comp, accum)
end

#Loading FullBinner
function load(
        ::Val{FullBinner{T, A}}, fn::Union{HDF5.File, HDF5.Group}, path::String) where {
        T, A <: AbstractArray}
    return FullBinner(collect.(eachcol(read_dataset(fn, "$(path)/values"))))
end

#Top-level load function
function load(fn::H5, path::String)
    observableType = eval(Meta.parse(read_dataset(fn, "$(path)/observableType")))
    return load(Val(observableType), fn, path)
end

"""
--------------------------------------------------------------------------------
Saving and loading Observables
--------------------------------------------------------------------------------
"""

"""
    function writeObservables!(fn::Union{HDF5.File,HDF5.Group}, obs::Observables)
Writes the observable type in a H5 object. 
This is done by calling save! on each field.
"""
function writeObservables!(fn::H5, obs::O) where {O <: AbstractObservables}
    o = create_group(fn, "observables")
    o["Type"] = String(Symbol(typeof(obs)))

    for field in fieldnames(O)
        currentObservableName = String(field)
        currentObservable = getfield(obs, field)
        save!(o, String(field), getfield(obs, field))
        o["$(currentObservableName)/observableType"] = String(Symbol(typeof(currentObservable)))
    end

    return nothing
end

"""
    function readObservables(fn::Union{HDF5.File,HDF5.Group})
Reads the observable type from a H5 object. 
It assumes that the appropriate constructor for the saved 
observable type exists. If no observable type is defined, 
it switches to the built-in observables.
"""
function readObservables(fn::H5)
    o = open_group(fn, "observables")
    observablesType = haskey(o, "Type") ? Symbol(read_dataset(o, "Type")) : nothing

    if observablesType === nothing
        return Observables([load(o, String(field)) for field in fieldnames(Observables)]...)
    else
        !(isdefined(Main, observablesType)) && error("$observablesType not defined")
        customObs = getfield(Main, observablesType)
        return customObs([load(o, String(field)) for field in fieldnames(customObs)]...)
    end
end

"""
--------------------------------------------------------------------------------
Saving and loading MonteCarlo
--------------------------------------------------------------------------------
"""

"""
    function writeMonteCarlo!(fn::H5, mc::MonteCarlo)
Writes the MonteCarlo object in a H5 object.
"""
function writeMonteCarlo!(fn::H5, mc::MonteCarlo)
    g = create_group(fn, "mc")
    writeLattice!(g, mc.lattice)
    writeMonteCarloParameters!(g, mc.parameters)
    writeMonteCarloStatistics!(g, mc.statistics)
    writeObservables!(g, mc.observables)
    return nothing
end

"""
    function writeMonteCarlo!(filename::String, mc::MonteCarlo)
Writes the MonteCarlo object into a HDF5 file.
"""
function writeMonteCarlo!(filename::String, mc::MonteCarlo)
    h5open(filename, "w") do f
        writeMonteCarlo!(f, mc)
    end
    return nothing
end

"""
    function readMonteCarlo(fn::H5)
Reads the MonteCarlo object from a H5 object.
"""
function readMonteCarlo(fn::H5)
    g = open_group(fn, "mc")
    lattice = readLattice(g)
    parameters = readMonteCarloParameters(g)
    statistics = readMonteCarloStatistics(g)
    observables = readObservables(g)
    return MonteCarlo(lattice, parameters, statistics, observables)
end

"""
    function readMonteCarlo(filename::String)
Reads the MonteCarlo object from a HDF5 file.
"""
function readMonteCarlo(filename::String)
    h5open(filename, "r") do f
        return readMonteCarlo(f)
    end
end

"""
--------------------------------------------------------------------------------
File read IO
--------------------------------------------------------------------------------
"""

function getFileNames(location::String; fileExtension = "h5")
    return filter(
        (x) -> contains(x, fileExtension) && !isdir(x), readdir(location, join = true))
end

function getFolderNames(location::String)
    return filter(isdir, readdir(location, join = true))
end

"""
--------------------------------------------------------------------------------
Functions for reading observable as measurements
--------------------------------------------------------------------------------
"""

function readObservable(f::H5, obs::Symbol)
    g = open_group(f, String(obs))
    return measurement(read_dataset(g, "mean"), read_dataset(g, "error"))
end

function readObservable(f::H5, obs::Vector{Symbol})
    observables = Vector{Measurement{Float64}}(undef, length(obs))
    for (i, o) in enumerate(obs)
        observables[i] = loadObservable(f, o)
    end
    return observables
end

function readObservable(fn::String, obs::Vector{Symbol})
    h5open(fn, "r") do f
        g = f["mc/observables"]
        return readObservable(g, obs)
    end
end

function readObservable(
        ::Val{:fixed_temperature}, fn::Vector{String}, obs::Vector{Symbol})
    observables = Matrix{Measurement}(undef, 1, length(obs))
    observables[1, :] .= mean(readObservable.(fn, Ref(obs)))
    return observables
end

function readObservable(
        ::Val{:changing_temperature}, fn::Vector{String}, obs::Vector{Symbol})
    observables = Matrix{Measurement}(undef, length(fn), length(obs))
    for (j, file) in enumerate(fn)
        observables[j, :] .= readObservable(file, obs)
    end
    return observables
end

function readObservable(fn::Vector{String}, obs::Vector{Symbol})
    β = map(i => h5read(i, "mc/parameters/beta"), fn)
    perm = sortperm(β)
    β = β[perm]
    fn = fn[perm]
    temperatures = unique(inv.(β))
    if isone(length(temperatures))
        return temperatures, readObservable(Val(:fixed_temperature), fn, obs)
    elseif length(temperatures) == length(fn)
        return temperatures, readObservable(Val(:changing_temperature), fn, obs)
    else
        error("Please limit files to (a) single β or (b) share no common β.")
    end
end

"""
--------------------------------------------------------------------------------
Process observables to store in a compact file format.
--------------------------------------------------------------------------------
"""

"""
function processObservables!(
        temperatures::Vector{Float64}, observables::Matrix{Measurement{T}},
        obs::Vector{Symbol}, f::H5) where {T}

Function to save an observables matrix into a H5 object.
"""
function processObservables!(
        temperatures::Vector{Float64}, observables::Matrix{Measurement{T}},
        obs::Vector{Symbol}, f::H5) where {T}
    f["temperatures"] = temperatures
    f["properties"] = String.(obs)
    g = create_group(f, "observables")
    for (i, o) in enumerate(obs)
        h = create_group(g, String(o))
        o = observables[:, i]
        h["val"] = getfield.(o, :val)
        h["err"] = getfield.(o, :err)
    end
end

"""
function processObservables!(fn::Vector{String}, obs::Vector{Symbol}, saveLocation::String)

Function that reads data from a vector of filenames and then 
saves them into a single HDF5 file. 
Can be used to collect a mean observable for a single β across runs 
or to collect an observable acrosss temperatures.
"""
function processObservables!(fn::Vector{String}, obs::Vector{Symbol}, saveLocation::String)
    temperatures, observables = readObservable(fn, obs)
    h5open(saveLocation, "w") do f
        processObservables!(temperatures, observables, obs, f)
    end
end

"""
function processObservables!(fn::Vector{Vector{String}}, obs::Vector{Symbol},
        saveLocation::String)
    
Function that reads data from a vector of vector of filenames and 
then saves them into a single HDF5 file.

Assumes that the files in the inner vector are a single run across temperatures 
and the outer vector is iterating over a collection of runs.

Please note that the collection of runs are assumed to have exactly the same 
β values.
"""
function processObservables!(fn::Vector{Vector{String}}, obs::Vector{Symbol},
        saveLocation::String)
    data = [readObservable(i, obs) for i in fn]
    temperatures = only(unique(getindex.(data, 1))) #Will error if the runs have different collection of β
    observables = getindex.(data, 2)
    meanObservables = Matrix{Measurement{Float64}}(undef, length(temperatures), length(obs))
    for i in 1:size(meanObservables, 1)
        for j in 1:size(meanObservables, 2)
            meanObservables[i, j] = mean([o[i, j] for o in observables])
        end
    end
    h5open(saveLocation, "w") do f
        processObservables!(temperatures, meanObservables, obs, f)
    end
end

"""
function loadProcessedObservables(fn::String)

Function that loads previously processed observables into memory. 
Please note that since this is a compact representation of data, 
all technical information about the run such as the amount of sweeps /
number of runs is no longer retained.
"""
function loadProcessedObservables(fn::String)
    h5open(fn) do f
        temperatures = read(f["temperatures"])
        properties = Symbol.(read(f["properties"]))
        g = open_group(f, "observables")
        observables = Matrix{Measurement{Float64}}(
            undef, length(temperatures), length(properties))
        for (i, o) in enumerate(properties)
            h = open_group(g, String(o))
            observables[:, i] = measurement.(read(h["val"]), read(h["err"]))
        end
        return temperatures, observables, properties
    end
end

"""
--------------------------------------------------------------------------------
Loading spin configurations
--------------------------------------------------------------------------------
"""

getSpins(fn::String) = Vector{SVector{3, Float64}}(eachcol(h5read(fn, "mc/lattice/spins")))

getNSites(fn::String) = size(h5read(fn, "mc/lattice/unitcell/basis"), 2)

function groupSpins(nSites, spins)
    groupedSpins = Vector{Vector{SVector{3, Float64}}}(undef, nSites)
    totalSpins = length(spins)
    for currentSite in 1:nSites
        groupedSpins[currentSite] = map(
            i -> getindex(spins, i), currentSite:nSites:totalSpins)
    end
    return groupedSpins
end