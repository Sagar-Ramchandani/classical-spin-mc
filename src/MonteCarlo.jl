"""
    const updateParameterDict
A dictionary used to lookup the initial parameters for 
spin update methods.
"""
const updateParameterDict = Dict(conicalUpdate => 1.0π)

"""
--------------------------------------------------------------------------------
Monte Carlo Structs
--------------------------------------------------------------------------------
"""

"""
    mutable struct MonteCarloStatistics
This is used to store information regarding the statistics gathered
during the Monte Carlo run.
"""
@kwdef mutable struct MonteCarloStatistics
    sweeps::Int = 0

    attemptedLocalUpdatesTotal::Int = 0
    acceptedLocalUpdatesTotal::Int = 0
    attemptedReplicaExchangesTotal::Int = 0
    acceptedReplicaExchangesTotal::Int = 0

    attemptedLocalUpdates::Int = 0
    acceptedLocalUpdates::Int = 0
    attemptedReplicaExchanges::Int = 0
    acceptedReplicaExchanges::Int = 0

    initializationTime::Float64 = time()
end

"""
    mutable struct MonteCarloParameters
This is used to store any parameters used to configure the Monte Carlo run.
"""
@kwdef mutable struct MonteCarloParameters{U <: AbstractRNG, UP <: Function}
    beta::Float64
    thermalizationSweeps::Int
    measurementSweeps::Int
    measurementRate::Int = 1
    microcanonicalRoundsPerSweep::Int = 0
    replicaExchangeRate::Int = 10
    randomizeInitialConfiguration::Bool = true
    reportInterval::Int = round(Int, 0.05 * (thermalizationSweeps + measurementSweeps))
    checkpointInterval::Int = 3600

    rng::U = copy(Random.GLOBAL_RNG)
    seed::UInt = rand(Random.RandomDevice(), UInt)
    sweep::Int = 0

    updateFunction::UP = sphericalUpdate
    updateParameter::Float64 = get(updateParameterDict, updateFunction, 0.0)
end

"""
    mutable struct MonteCarlo
This is used to store all information about a Monte Carlo run.
"""
mutable struct MonteCarlo{T <: Lattice, P <: MonteCarloParameters, O <: AbstractObservables}
    lattice::T
    parameters::P
    statistics::MonteCarloStatistics
    observables::O

    function MonteCarlo(lattice::T, parameters::P, statistics::MonteCarloStatistics,
            observables::O) where {
            T <: Lattice, P <: MonteCarloParameters, O <: AbstractObservables}
        mc = new{T, P, O}(deepcopy(lattice), parameters, statistics, observables)
        Random.seed!(mc.parameters.rng, mc.parameters.seed)
        return mc
    end
end

"""
    mutable struct MonteCarloAnnealing
This is a wrapper around MonteCarlo and is used to store the 
multiple Monte Carlo structs created during a simulated annealing run.
"""
mutable struct MonteCarloAnnealing{}
    MonteCarloObjects::Vector{MonteCarlo}
end

"""
    mutable struct MonteCarloExchange
This is a wrapper around MonteCarlo and is used to store the 
multiple Monte Carlo structs and resources needed for distributed 
computing during a Replica Exchange (Parallel tempering) run.
"""
mutable struct MonteCarloExchange{T, M, A <: AbstractVector{MonteCarlo}}
    MonteCarloObjects::A
    betas::Vector{Float64}
    channelsUp::Vector{RemoteChannel{Channel{T}}}
    channelsDown::Vector{RemoteChannel{Channel{T}}}
    exchangeMethod::Val{M}
end

"""
--------------------------------------------------------------------------------
Monte Carlo Constructors
--------------------------------------------------------------------------------
"""

"""
    function MonteCarloParameters(
    beta::Float64,
    thermalizationSweeps::Int,
    measurementSweeps::Int)

Constructor that uses reasonable defaults for MonteCarloParameters.
"""
function MonteCarloParameters(
        beta::Float64,
        thermalizationSweeps::Int,
        measurementSweeps::Int)
    return MonteCarloParameters(beta = beta, thermalizationSweeps = thermalizationSweeps,
        measurementSweeps = measurementSweeps)
end

"""
    function MonteCarlo(lattice::T, parameters::MonteCarloParameters,
    observables::Observables) where {T<:Lattice}

Constructor that generates a MonteCarloStatistics struct automatically.
"""
function MonteCarlo(lattice::T, parameters::MonteCarloParameters,
        observables::O) where {T <: Lattice, O <: AbstractObservables}
    return MonteCarlo(lattice, parameters, MonteCarloStatistics(), observables)
end

"""
    function MonteCarlo(lattice::T, parameters::MonteCarloParameters) where {T<:Lattice}

Constructor that generates a MonteCarloStatistics struct and built-in 
Observables struct automatically.
"""
function MonteCarlo(lattice::T, parameters::MonteCarloParameters) where {T <: Lattice}
    return MonteCarlo(lattice, parameters, MonteCarloStatistics(), Observables(lattice))
end

"""
    function MonteCarloAnnealing(mc::MonteCarlo, betas::Vector{Float64})

Constructor for generating a MonteCarloAnnealing struct used for 
Simulated Annealing runs.
"""
function MonteCarloAnnealing(mc::MonteCarlo, betas::Vector{Float64})
    simulations = Vector{MonteCarlo}(undef, length(betas))

    #Check if betas are sorted in ascending order else sort them.
    if !(issorted(betas))
        @warn "Input βs are not sorted. Sorting them anyways."
        sort!(betas)
    end

    for (i, beta) in enumerate(betas)
        #create one simulation for each provided beta based on the specified mc template
        simulations[i] = deepcopy(mc)
        simulations[i].parameters.beta = beta
        if i != 1
            #if this is not the first simulation, set spin randomization false
            simulations[i].parameters.randomizeInitialConfiguration = false
        end
    end
    return MonteCarloAnnealing(simulations)
end

"""
    function MonteCarloExchange(mc::MonteCarlo, betas::Vector{Float64})

Constructor for generating a MonteCarloExchange struct used for 
Replica Exchange (Parallel tempering) runs.

Available methods are :parallel and :serial
:serial offers better stability at the cost of performance.


!!! warning "Available workers" 
    At the present time, this module uses workers 1:length(betas).
    Please verify that these workers are created and available before 
    creating this struct. At a future time, arbritary workers and 
    their automatic creation may be supported.
"""

function MonteCarloExchange(mc::MonteCarlo, betas::Vector{Float64}; method = :parallel)
    simulations = Vector{MonteCarlo}(undef, length(betas))

    #Check if betas are sorted in ascending order else sorted them
    if !(issorted(betas))
        @warn "Input βs are not sorted. Sorting them anyways."
        sort!(betas)
    end

    for (i, beta) in enumerate(betas)
        #create one simulation for each provided beta based on the specified mc template
        simulations[i] = deepcopy(mc)
        simulations[i].parameters.beta = beta
    end
    return MonteCarloExchange(simulations, betas, createChannels()..., Val(method))
end

"""
--------------------------------------------------------------------------------
Extend Base for Monte Carlo structs
--------------------------------------------------------------------------------
"""

function Base.:show(io::IO, mc::MonteCarlo{T}) where {T <: Lattice}
    println(io,
        "MonteCarlo with β=$(mc.parameters.beta) and update function: $(mc.parameters.updateFunction)")
end

function Base.:show(io::IO, parameters::MonteCarloParameters)
    println(io,
        "MonteCarloParameters with β=$(parameters.beta) and update function: $(parameters.updateFunction)")
end

function Base.:show(io::IO, statistics::MonteCarloStatistics)
    time = Dates.format(unix2datetime(statistics.initializationTime), "dd u yyyy HH:MM:SS")
    str = "MonteCarloStatistics with $(statistics.sweeps) Sweeps, initialized at $time\n"
    if !(statistics.attemptedLocalUpdatesTotal == 0)
        str *= @sprintf("\tUpdate acceptance rate: %.2f%%\n",
            100 *
            statistics.acceptedLocalUpdatesTotal/statistics.attemptedLocalUpdatesTotal)
    end
    if !(statistics.attemptedReplicaExchangesTotal == 0)
        str *= @sprintf("\tReplica acceptance rate: %.2f%%\n",
            100 *
            statistics.acceptedReplicaExchangesTotal/statistics.attemptedReplicaExchangesTotal)
    end
    println(
        io,
        str
    )
end

import Base: ==

function ==(a::T, b::T) where {T <: Union{MonteCarloStatistics, MonteCarloParameters}}
    fields = fieldnames(T)
    status = true
    for field in fields
        status = status && (getfield(a, field) == getfield(b, field))
    end
    return status
end

function ==(a::T, b::T) where {T <: MonteCarlo}
    fields = fieldnames(T)
    status = true
    for field in fields
        status = status && (getfield(a, field) == getfield(b, field))
        if status == false
            println(field)
        end
    end
    return status
end

"""
--------------------------------------------------------------------------------
initSpinConfiguration functions
--------------------------------------------------------------------------------
"""

"""
    function initSpinConfiguration!(lattice::Lattice{D,N}, f::typeof(conicalUpdate), rng=Random.GLOBAL_RNG) where {D,N}
Initialize the spin configuration of a lattice using the conicalUpdate method.

!!! warning "Pre-initialize for conical updates"
    Since conicalUpdate forms a cone around a point on a unit sphere, it
    cannot be used to reliably initialize a spin configuration unless it was already initialized 
    using a different method first. This method is provided only for consistency.
"""
function initSpinConfiguration!(lattice::Lattice{D, N}, f::typeof(conicalUpdate),
        rng = Random.GLOBAL_RNG) where {D, N}
    @warn "Conical updates only work if the lattice is already initialized,
    using sphericalUpdate instead"
    for i in 1:length(lattice)
        setSpin!(lattice, i, sphericalUpdate(rng))
    end
end

"""
    function initSpinConfiguration!(lattice::Lattice{D,N}, f::Function, rng=Random.GLOBAL_RNG) where {D,N}
Initialize the spin configuration of a lattice using the passed function.
"""
function initSpinConfiguration!(
        lattice::Lattice{D, N}, f::Function, rng = Random.GLOBAL_RNG) where {D, N}
    for i in 1:length(lattice)
        setSpin!(lattice, i, f(rng))
    end
end

"""
    function initSpinConfiguration!(mc::MonteCarlo{T}) where {T<:Lattice}
Initialize the spin configuration for a Monte Carlo struct by using the 
function specified in MonteCarloParameters.
"""
function initSpinConfiguration!(mc::MonteCarlo{T}) where {T <: Lattice}
    if (mc.parameters.sweep == 0) && mc.parameters.randomizeInitialConfiguration
        initSpinConfiguration!(mc.lattice, mc.parameters.updateFunction, mc.parameters.rng)
        mc.parameters.randomizeInitialConfiguration = false
    end
end

"""
--------------------------------------------------------------------------------
localSweep functions
--------------------------------------------------------------------------------
"""

"""
    function localUpdate(mc::MonteCarlo{T,P}, proposalSite::Int64, newSpinState::SVector{3,Float64}) where {T<:Lattice,P<:MonteCarloParameters}
Performs a local update on the current spin configuration, using the proposed new spin state.
"""
function localUpdate(mc::MonteCarlo{T, P}, proposalSite::Int64,
        newSpinState::SVector{3, Float64}) where {T <: Lattice, P <: MonteCarloParameters}
    energyDifference = getEnergyDifference(mc.lattice, proposalSite, newSpinState)
    #check acceptance of new configuration
    mc.statistics.attemptedLocalUpdates += 1
    p = exp(-mc.parameters.beta * energyDifference)
    if (rand(mc.parameters.rng) < min(1.0, p))
        setSpin!(mc.lattice, proposalSite, newSpinState)
        mc.statistics.acceptedLocalUpdates += 1
    else
        energyDifference = 0.0
    end
    return energyDifference
end

#Note: Possibly change this function to use the passed function instead of referring to MCP.
"""
    function localSweep(::Function, mc::MonteCarlo{T}, energy::Float64) where {T<:Lattice}
Performs a local sweep on the current spin configuration, using the function specified in 
MonteCarloParameters. This covers only the case where the update function
takes the rng as it's only parameter.
"""
function localSweep(::Function, mc::MonteCarlo{T}, energy::Float64) where {T <: Lattice}
    for _ in 1:length(mc.lattice)
        site = rand(mc.parameters.rng, 1:length(mc.lattice))
        newSpinState = mc.parameters.updateFunction(mc.parameters.rng)
        energy += localUpdate(mc, site, newSpinState)
    end
    return energy
end

"""
    function localSweep(::typeof(conicalUpdate), mc::MonteCarlo{T}, energy::Float64) where {T<:Lattice}
Performs a local sweep on the current spin configuration using the conicalUpdate method. 
This is a specialization of the method and not the general case.
"""
function localSweep(
        ::typeof(conicalUpdate), mc::MonteCarlo{T}, energy::Float64) where {T <: Lattice}
    for _ in 1:length(mc.lattice)
        site = rand(mc.parameters.rng, 1:length(mc.lattice))
        newSpinState = mc.parameters.updateFunction(
            getSpin(mc.lattice, site), mc.parameters.updateParameter, mc.parameters.rng)
        energy += localUpdate(mc, site, newSpinState)
    end

    #Conical update adapt parameter (cone angle)
    adaptiveFactor = 0.5 / (1 - (mc.statistics.acceptedLocalUpdates /
                       mc.statistics.attemptedLocalUpdates))
    mc.parameters.updateParameter *= adaptiveFactor
    if mc.parameters.updateParameter > 1.0π
        mc.parameters.updateParameter = 1.0π
    end
    return energy
end

"""
    function localSweep(mc::MonteCarlo{T}, energy::Float64) where {T<:Lattice}
Performs a local sweep on the current spin configuration using the method specified
in MonteCarloParameters. This is the top-level method.
"""
function localSweep(mc::MonteCarlo{T}, energy::Float64) where {T <: Lattice}
    energy = localSweep(mc.parameters.updateFunction, mc, energy)
    mc.parameters.sweep += 1
    mc.statistics.sweeps += 1
    return energy
end

"""
--------------------------------------------------------------------------------
microcanonicalSweep functions
--------------------------------------------------------------------------------
"""

"""
    function microcanonicalSweep!(lattice::Lattice{D,N}, rounds::Int, rng=Random.GLOBAL_RNG) where {D,N}
Performs a microcanonical (overrelaxation) sweep on the current spin configuration.
This is particularly useful for low-temperature simulations where the spins get stuck.
"""
function microcanonicalSweep!(
        lattice::Lattice{D, N}, rounds::Int, rng = Random.GLOBAL_RNG) where {D, N}
    basisLength = length(lattice.unitcell)
    for _ in 1:rounds
        #Deterministic first iteration
        for i in 1:basisLength
            sublatticeOrdered = range(i, length(lattice), step = basisLength)
            for site in sublatticeOrdered
                newSpinState = microcanonicalRotation(lattice, site)
                setSpin!(lattice, site, newSpinState)
            end
        end
        #Random second iteration
        for i in 1:basisLength
            sublatticeOrdered = range(i, length(lattice), step = basisLength)
            for site in sublatticeOrdered
                setSpin!(lattice, site, microcanonicalRotationRandom(lattice, site, rng))
            end
        end
    end
    return nothing
end

"""
    function microcanonicalSweep!(mc::MonteCarlo{T}) where {T<:Lattice}
Performs a microcanonical (overrelaxation) sweep on current spin configuration,
based on parameters from MonteCarloParameters.
"""
function microcanonicalSweep!(mc::MonteCarlo{T}) where {T <: Lattice}
    return microcanonicalSweep!(
        mc.lattice, mc.parameters.microcanonicalRoundsPerSweep, mc.parameters.rng)
end

"""
--------------------------------------------------------------------------------
replicaExchange functions
--------------------------------------------------------------------------------
"""

"""
    function replicaExchange!(mc::MonteCarlo{T}, energy::Float64, betas::Vector{Float64}, channelsUp, channelsDown, label) where {T<:Lattice}
Performs a parallel replica exchange sweep based on parameters from the MonteCarlo object on the current worker.
"""
function replicaExchange!(
        ::Val{:parallel}, mc::MonteCarlo{T}, energy::Float64, betas::Vector{Float64},
        channelsUp, channelsDown, label) where {T <: Lattice}
    #determine replica partner rank
    rank = myid() - 1
    numberWorkers = length(betas)
    if iseven(mc.parameters.sweep ÷ mc.parameters.replicaExchangeRate)
        partnerRank = iseven(rank) ? rank + 1 : rank - 1
    else
        partnerRank = iseven(rank) ? rank - 1 : rank + 1
    end
    if partnerRank >= 1 && partnerRank <= numberWorkers
        if partnerRank > rank
            chPut = channelsUp[2]
            chTake = channelsDown[2]
        else
            chPut = channelsDown[1]
            chTake = channelsUp[1]
        end
        #put own energy on remoteChannel
        put!(chPut, energy)
        #Take partner energy from remoteChannel
        partnerEnergy = take!(chTake)

        #check acceptance of new configuration
        mc.statistics.attemptedReplicaExchanges += 1
        exchangeAccepted = false

        if iseven(rank)
            p = exp(-(betas[rank] - betas[partnerRank]) * (partnerEnergy - energy))
            exchangeAccepted = (rand(mc.parameters.rng) < min(1.0, p)) ? true : false
            put!(chPut, exchangeAccepted)
        else
            exchangeAccepted = take!(chTake)
        end
        if (exchangeAccepted)
            energy = partnerEnergy
            put!(chPut, mc.lattice.spins)
            mc.lattice.spins = take!(chTake)
            put!(chPut, label)
            label = take!(chTake)
            if rank == numberWorkers
                label = -1
            end
            if rank == 1
                label = 1
            end

            mc.statistics.acceptedReplicaExchanges += 1
        end
    end
    return (energy, label)
end

"""
    function replicaExchange!(mc::MonteCarlo{T}, energy::Float64, betas::Vector{Float64}, channelsUp, channelsDown, label) where {T<:Lattice}
Performs a serial replica exchange sweep based on parameters from the MonteCarlo object on the current worker.
"""
function replicaExchange!(
        ::Val{:serial}, mc::MonteCarlo{T}, energy::Float64, betas::Vector{Float64},
        channelsUp, channelsDown, label) where {T <: Lattice}
    #determine replica partner rank
    rank = myid() - 1
    numberWorkers = length(betas)
    if rank == numberWorkers
        mc.statistics.attemptedReplicaExchanges += 1
        partnerRank = rank - 1
        chPut = channelsDown[1]
        chTake = channelsUp[1]
        partnerEnergy = take!(chTake)
        p = exp((betas[partnerRank] - betas[rank]) * (partnerEnergy - energy))
        exchangeAccepted = (rand(mc.parameters.rng) < min(1.0, p)) ? true : false
        put!(chPut, exchangeAccepted)
        if exchangeAccepted
            put!(chPut, energy)
            put!(chPut, mc.lattice.spins)
            energy = partnerEnergy
            mc.lattice.spins = take!(chTake)
            mc.statistics.acceptedReplicaExchanges += 1
        end
    else
        mc.statistics.attemptedReplicaExchanges += 1
        partnerRank = rank + 1
        chPut = channelsUp[2]
        chTake = channelsDown[2]
        put!(chPut, energy)
        exchangeAccepted = take!(chTake)
        if exchangeAccepted
            put!(chPut, mc.lattice.spins)
            energy = take!(chTake)
            mc.lattice.spins = take!(chTake)
            mc.statistics.acceptedReplicaExchanges += 1
        end
        if !(rank == 1)
            mc.statistics.attemptedReplicaExchanges += 1
            partnerRank = rank - 1
            chPut = channelsDown[1]
            chTake = channelsUp[1]
            partnerEnergy = take!(chTake)
            p = exp(-(betas[rank] - betas[partnerRank]) * (partnerEnergy - energy))
            exchangeAccepted = (rand(mc.parameters.rng) < min(1.0, p)) ? true : false
            put!(chPut, exchangeAccepted)
            if exchangeAccepted
                put!(chPut, energy)
                put!(chPut, mc.lattice.spins)
                energy = partnerEnergy
                mc.lattice.spins = take!(chTake)
                mc.statistics.acceptedReplicaExchanges += 1
            end
        end
    end
    return (energy, 0)
end

"""
--------------------------------------------------------------------------------
Monte Carlo Functions
--------------------------------------------------------------------------------
"""

"""
    function updateTotalStatistics!(statistics::MonteCarloStatistics)
Updates the running total of statistics while resetting statistics
from the current sweep.
"""
function updateTotalStatistics!(statistics::MonteCarloStatistics)
    #update running total of statistics
    statistics.attemptedLocalUpdatesTotal += statistics.attemptedLocalUpdates
    statistics.acceptedLocalUpdatesTotal += statistics.acceptedLocalUpdates
    statistics.attemptedReplicaExchangesTotal += statistics.attemptedReplicaExchanges
    statistics.acceptedReplicaExchangesTotal += statistics.acceptedReplicaExchanges

    #Reset current statistics
    statistics.attemptedLocalUpdates = 0
    statistics.acceptedLocalUpdates = 0
    statistics.attemptedReplicaExchanges = 0
    statistics.acceptedReplicaExchanges = 0
    return nothing
end

"""
    function createChannels()
Creates the RemoteChannels required for worker-worker communication used in 
Replica Exchange (Parallel Tempering) runs.
"""
function createChannels()
    channelsUp = [fetch(@spawnat i RemoteChannel(() -> Channel{Any}(1), myid()))
                  for i in procs()]
    channelsDown = circshift(
        [fetch(@spawnat i RemoteChannel(() -> Channel{Any}(1), myid())) for i in procs()],
        -1)
    return (channelsUp, channelsDown)
end

"""
    function printStatistics!(mc::MonteCarlo{T}; replica=false) where {T<:Lattice}
Print statistics from the current status of the MonteCarlo run.
"""
function printStatistics!(mc::MonteCarlo{T}; replica = false) where {T <: Lattice}
    t = time()
    if mc.parameters.sweep % mc.parameters.reportInterval == 0
        #collect statistics
        totalSweeps = mc.parameters.thermalizationSweeps + mc.parameters.measurementSweeps
        progress = 100.0 * mc.parameters.sweep / totalSweeps
        thermalized = (mc.parameters.sweep >= mc.parameters.thermalizationSweeps) ? "YES" :
                      "NO"
        sweeprate = mc.statistics.sweeps / (t - mc.statistics.initializationTime)
        sweeptime = 1.0 / sweeprate
        eta = (totalSweeps - mc.parameters.sweep) / sweeprate
        localUpdateAcceptanceRate = 100.0 * mc.statistics.acceptedLocalUpdates /
                                    mc.statistics.attemptedLocalUpdates

        #print statistics
        str = ""
        str *= @sprintf("Sweep %d / %d (%.1f%%)", mc.parameters.sweep, totalSweeps,
            progress)
        str *= @sprintf("\t\tETA : %s\n",
            Dates.format(
                Dates.now() + Dates.Second(round(Int64, eta)), "dd u yyyy HH:MM:SS"))
        str *= @sprintf("\t\tthermalized : %s\n", thermalized)
        str *= @sprintf("\t\tsweep rate : %.1f sweeps/s\n", sweeprate)
        str *= @sprintf("\t\tsweep duration : %.3f ms\n", sweeptime*1000)
        str *= @sprintf("\t\tupdate acceptance rate: %.2f%%\n", localUpdateAcceptanceRate)
        str *= @sprintf("\t\tupdate parameter: %.3f \n", mc.parameters.updateParameter)
        if replica
            str *= @sprintf("\t\treplica acceptanced : %.3f \n",
                mc.statistics.acceptedReplicaExchanges)
            str *= @sprintf("\t\treplica attempted : %.3f \n",
                mc.statistics.attemptedReplicaExchanges)
        end
        str *= @sprintf("\n")
        print(str)

        #update total statistics and reset current statistics
        updateTotalStatistics!(mc.statistics)
    end
end

"""
    function sanityChecks(mc::MonteCarlo{T}, outfile::Union{String,Nothing}=nothing) where {T<:Lattice}
Perform sanity checks on a MonteCarlo struct to prevent undefined behaviour during the run.
"""
function sanityChecks(
        mc::MonteCarlo{T}, outfile::Union{String, Nothing} = nothing) where {T <: Lattice}
    #init IO
    enableOutput = typeof(outfile) != Nothing
    if enableOutput
        isfile(outfile) && error("File ", outfile, " already exists. Terminating.")
    end

    #Check if microcanonical conditions apply
    if mc.parameters.microcanonicalRoundsPerSweep != 0
        for site in 1:length(mc.lattice)
            if getInteractionOnsite(mc.lattice, site) != @SMatrix zeros(3, 3)
                error("Microcanonical updates are only supported for models without on-size interactions.")
            end
        end
    end

    #Update initializationTime
    mc.statistics.initializationTime = time()
    return enableOutput
end

"""
--------------------------------------------------------------------------------
Monte Carlo run! Functions
--------------------------------------------------------------------------------
"""

"""
    function run!(mc::MonteCarlo{T}; outfile::Union{String,Nothing}=nothing) where {T<:Lattice}
Dispatches run to perform a fixed-temperature Monte Carlo run.
"""
function run!(
        mc::MonteCarlo{T}; outfile::Union{String, Nothing} = nothing) where {T <: Lattice}
    enableOutput = sanityChecks(mc, outfile)

    #init spin configuration
    initSpinConfiguration!(mc)

    #init Monte Carlo run
    totalSweeps = mc.parameters.thermalizationSweeps + mc.parameters.measurementSweeps
    energy = getEnergy(mc.lattice)

    #launch Monte Carlo run
    lastCheckpointTime = time()
    @printf("Simulation started on %s.\n\n",
        Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))

    while mc.parameters.sweep < totalSweeps
        #perform local sweep
        energy = localSweep(mc, energy)
        #perform microcanonical sweep
        if (mc.parameters.microcanonicalRoundsPerSweep != 0) &&
           (mc.parameters.sweep % mc.parameters.microcanonicalRoundsPerSweep) == 0
            microcanonicalSweep!(mc)
        end
        #perform measurement
        if mc.parameters.sweep >= mc.parameters.thermalizationSweeps
            if mc.parameters.sweep % mc.parameters.measurementRate == 0
                performMeasurements!(mc.observables, mc.lattice, energy)
            end
        end

        #runtime statistics
        printStatistics!(mc)

        #write checkpoint
        if enableOutput
            checkpointPending = time() - lastCheckpointTime >=
                                mc.parameters.checkpointInterval
            if checkpointPending
                writeCheckpoint!(outfile, mc)
                lastCheckpointTime = time()
                @printf("Checkpoint written on %s.\n",
                    Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))
            end
        end
    end

    #perform post-measurement
    performPostMeasurements!(mc.observables, mc.lattice, mc.parameters.beta)

    #write final checkpoint
    if enableOutput
        writeMonteCarlo!(outfile, mc)
        @printf("Checkpoint written on %s.\n",
            Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))
    end

    #return
    @printf("Simulation finished on %s.\n", Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))
    return nothing
end

"""
    function anneal!(mc::M, betas::Vector{Float64}) where {M<:MonteCarlo}
Anneals a given MonteCarlo object through the given array of betas 
without creating intermediate MonteCarloAnnealing and
MonteCarlo objects.
"""
function anneal!(mc::M, betas::Vector{Float64}) where {M <: MonteCarlo}
    measurementSweeps = mc.parameters.measurementSweeps
    measurementSweeps = 0
    for (i, b) in betas
        mc.parameters.sweeps = 0
        mc.parameters.beta = b
        if i != 1
            mc.parameters.randomizeInitialConfiguration = false
        end
        if i == n
            mc.parameters.measurementSweeps = measurementSweeps
        end
        run!(mc, outfile = nothing)
    end
end

"""
    function run!(mcs::MonteCarloAnnealing; outfile::Union{String,Nothing}=nothing)
Dispatches run to perform a Simulated Annealing run.
"""
function run!(mcs::MonteCarloAnnealing; outfile::Union{String, Nothing} = nothing)
    for (i, mc) in enumerate(mcs.MonteCarloObjects)
        if i != 1
            #If not the first simulation, copy spin configuration from the previous one
            mcs.MonteCarloObjects[i].lattice = deepcopy(mcs.MonteCarloObjects[i - 1]).lattice
        end
        out = outfile === nothing ? outfile : outfile * "." * string(i - 1)
        run!(mc, outfile = out)
    end
end

"""
    function run!(mcs::MonteCarloExchange, outfile::Union{String,Nothing}=nothing)
Dispatches run to perform a Replica Exchange (Parallel Tempering) run.
"""
function run!(mcs::MonteCarloExchange, outfile::Union{String, Nothing} = nothing)
    workersList = workers()
    function sim(i)
        return @spawnat workersList[i] run!(
            mcs.MonteCarloObjects[i], mcs.betas, mcs.exchangeMethod,
            mcs.channelsUp[i:(i + 1)], mcs.channelsDown[i:(i + 1)],
            (outfile === nothing ? outfile : outfile * "." * string(i - 1)))
    end
    results = map(sim, eachindex(workersList))
    mcs.MonteCarloObjects = fetch.(results)
end

"""
    function run!(mc::MonteCarlo{T}, betas::Vector{Float64}, channelsUp::Vector{RemoteChannel{Channel{C}}},
    channelsDown::Vector{RemoteChannel{Channel{C}}},
Internal function used for parallel tempering runs.
"""
function run!(mc::MonteCarlo{T}, betas::Vector{Float64}, method,
        channelsUp::Vector{RemoteChannel{Channel{C}}},
        channelsDown::Vector{RemoteChannel{Channel{C}}},
        outfile::Union{String, Nothing} = nothing) where {T <: Lattice, C}
    enableOutput = sanityChecks(mc, outfile)
    rank = myid() - 1
    nSimulations = nworkers()
    labels = Vector{Int}(undef,
        floor(Int,
            (mc.parameters.thermalizationSweeps + mc.parameters.measurementSweeps) /
            mc.parameters.replicaExchangeRate) + 1)
    if rank == 1
        currentLabel = 1
    elseif rank == nSimulations
        currentLabel = -1
    else
        currentLabel = 0
    end
    labels[1] = currentLabel

    println("Running replica exchanges across $(nSimulations) simulations")

    #init spin configuration
    initSpinConfiguration!(mc)

    #init Monte Carlo run
    totalSweeps = mc.parameters.thermalizationSweeps + mc.parameters.measurementSweeps
    energy = getEnergy(mc.lattice)

    #launch Monte Carlo run
    lastCheckpointTime = time()
    if rank == 1
        @printf("Simulation started on %s.\n\n",
            Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))
    end

    while mc.parameters.sweep < totalSweeps
        #perform local sweep
        energy = localSweep(mc, energy)

        #perform microcanonical sweep
        if (mc.parameters.microcanonicalRoundsPerSweep != 0) &&
           (mc.parameters.sweep % mc.parameters.microcanonicalRoundsPerSweep) == 0
            microcanonicalSweep!(mc)
        end

        #perform replica exchange
        if mc.parameters.sweep % mc.parameters.replicaExchangeRate == 0
            energy, currentLabel = replicaExchange!(method,
                mc, energy, betas, channelsUp, channelsDown, currentLabel)
            labels[floor(Int, mc.parameters.sweep / mc.parameters.replicaExchangeRate) + 1] = currentLabel
            push!(mc.observables.labels, currentLabel)
        end

        #perform measurement
        if mc.parameters.sweep >= mc.parameters.thermalizationSweeps
            if mc.parameters.sweep % mc.parameters.measurementRate == 0
                performMeasurements!(mc.observables, mc.lattice, energy)
            end
        end

        #runtime statistics
        if rank == 1
            printStatistics!(mc, replica = true)
        end

        #write checkpoint
        if enableOutput
            checkpointPending = time() - lastCheckpointTime >=
                                mc.parameters.checkpointInterval
            if checkpointPending
                writeCheckpoint!(outfile, mc)
                lastCheckpointTime = time()
                rank == 1 && @printf("Checkpoint written on %s.\n",
                    Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))
            end
        end
    end

    #perform post-measurement
    performPostMeasurements!(mc.observables, mc.lattice, mc.parameters.beta)

    #write final checkpoint
    if enableOutput
        writeMonteCarlo!(outfile, mc)
        rank == 1 && @printf("Checkpoint written on %s.\n",
            Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))
    end

    #return
    rank == 1 && @printf("Simulation finished on %s.\n",
        Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))
    return mc
end
