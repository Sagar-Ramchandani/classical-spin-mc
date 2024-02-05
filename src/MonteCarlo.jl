using Random
using Dates
using Printf
using Distributed

const updateParameterDict = Dict(conicalUpdate => 1.0π)

"""
--------------------------------------------------------------------------------
Monte Carlo Structs
--------------------------------------------------------------------------------
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

@kwdef mutable struct MonteCarloParameters{U<:AbstractRNG,UP<:Function}
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

mutable struct MonteCarlo{T<:Lattice,P<:MonteCarloParameters}
    lattice::T
    parameters::P
    statistics::MonteCarloStatistics
    observables::Observables

    function MonteCarlo(lattice::T, parameters::P, statistics::MonteCarloStatistics,
        observables::Observables) where {T<:Lattice,P<:MonteCarloParameters}
        mc = new{T,P}(deepcopy(lattice), parameters, statistics, observables)
        Random.seed!(mc.parameters.rng, mc.parameters.seed)
        return mc
    end
end

mutable struct MonteCarloAnnealing{}
    MonteCarloObjects::Vector{MonteCarlo}
end

mutable struct MonteCarloExchange{T}
    MonteCarloObjects::Vector{MonteCarlo}
    betas::Vector{Float64}
    channelsUp::Vector{RemoteChannel{Channel{T}}}
    channelsDown::Vector{RemoteChannel{Channel{T}}}
end

"""
--------------------------------------------------------------------------------
Monte Carlo Constructors
--------------------------------------------------------------------------------
"""

function MonteCarloParameters(
    beta::Float64,
    thermalizationSweeps::Int,
    measurementSweeps::Int)

    return MonteCarloParameters(beta=beta, thermalizationSweeps=thermalizationSweeps, measurementSweeps=measurementSweeps)
end


function MonteCarlo(lattice::T, parameters::MonteCarloParameters,
    observables::Observables) where {T<:Lattice}
    return MonteCarlo(lattice, parameters, MonteCarloStatistics(), observables)
end

function MonteCarlo(lattice::T, parameters::MonteCarloParameters) where {T<:Lattice}
    return MonteCarlo(lattice, parameters, MonteCarloStatistics(), Observables(lattice))
end

function MonteCarloAnnealing(mc::MonteCarlo, betas::Vector{Float64})
    simulations = Vector{MonteCarlo}(undef, length(betas))

    #Check if betas are sorted in descending order else sorted them
    if !(issorted(betas, rev=true))
        @warn "Input βs are not sorted. Sorting them anyways."
        sort!(betas, rev=true)
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

function MonteCarloExchange(mc::MonteCarlo, betas::Vector{Float64})
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
    return MonteCarloExchange(simulations, betas, createChannels()...)
end

"""
--------------------------------------------------------------------------------
Extend Base for Monte Carlo structs
--------------------------------------------------------------------------------
"""

function Base.:show(io::IO, mc::MonteCarlo{T}) where {T<:Lattice}
    println(io, "MonteCarlo with β=$(mc.parameters.beta) and update function: $(mc.parameters.updateFunction)")
end

function Base.:show(io::IO, parameters::MonteCarloParameters)
    println(io, "MonteCarloParameters with β=$(parameters.beta) and update function: $(parameters.updateFunction)")
end

function Base.:show(io::IO, statistics::MonteCarloStatistics)
    time = Dates.format(unix2datetime(statistics.initializationTime), "dd u yyyy HH:MM:SS")
    str = "MonteCarloStatistics with $(statistics.sweeps) Sweeps, initialized at $time\n"
    if !(statistics.attemptedLocalUpdatesTotal == 0)
        str *= @sprintf("\tUpdate acceptance rate: %.2f%%\n", 100 * statistics.acceptedLocalUpdatesTotal / statistics.attemptedLocalUpdatesTotal)
    end
    if !(statistics.attemptedReplicaExchangesTotal == 0)
        str *= @sprintf("\tReplica acceptance rate: %.2f%%\n", 100 * statistics.acceptedReplicaExchangesTotal / statistics.attemptedReplicaExchangesTotal)
    end
    println(
        io,
        str
    )
end

import Base: ==

function ==(a::T, b::T) where {T<:Union{MonteCarloStatistics,MonteCarloParameters}}
    fields = fieldnames(T)
    status = true
    for field in fields
        status = status && (getfield(a, field) == getfield(b, field))
    end
    return status
end

function ==(a::T, b::T) where {T<:MonteCarlo}
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

function initSpinConfiguration!(lattice::Lattice{D,N}, f::typeof(conicalUpdate), rng=Random.GLOBAL_RNG) where {D,N}
    @warn "Conical updates only work if the lattice is already initialized"
    for i in 1:length(lattice)
        setSpin!(lattice, i, f(getSpin(lattice, i), 1.0π, rng))
    end
end

function initSpinConfiguration!(lattice::Lattice{D,N}, f::Function, rng=Random.GLOBAL_RNG) where {D,N}
    for i in 1:length(lattice)
        setSpin!(lattice, i, f(rng))
    end
end

function initSpinConfiguration!(mc::MonteCarlo{T}) where {T<:Lattice}
    if (mc.parameters.sweep == 0) && mc.parameters.randomizeInitialConfiguration
        initSpinConfiguration!(mc.lattice, mc.parameters.updateFunction, mc.parameters.rng)
    end
end

"""
--------------------------------------------------------------------------------
localSweep functions
--------------------------------------------------------------------------------
"""

function localUpdate(mc::MonteCarlo{T,P}, proposalSite::Int64, newSpinState::SVector{3,Float64}) where {T<:Lattice,P<:MonteCarloParameters}
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

function localSweep(::Function, mc::MonteCarlo{T}, energy::Float64) where {T<:Lattice}
    for _ in 1:length(mc.lattice)
        site = rand(mc.parameters.rng, 1:length(mc.lattice))
        newSpinState = mc.parameters.updateFunction(mc.parameters.rng)
        energy += localUpdate(mc, site, newSpinState)
    end
    return energy
end

function localSweep(::typeof(conicalUpdate), mc::MonteCarlo{T}, energy::Float64) where {T<:Lattice}
    for _ in 1:length(mc.lattice)
        site = rand(mc.parameters.rng, 1:length(mc.lattice))
        newSpinState = mc.parameters.updateFunction(getSpin(mc.lattice, site), mc.parameters.updateParameter, mc.parameters.rng)
        energy += localUpdate(mc, site, newSpinState)
    end

    #Conical update adapt parameter (cone angle)
    adaptiveFactor = 0.5 / (1 - (mc.statistics.acceptedLocalUpdates / mc.statistics.attemptedLocalUpdates))
    mc.parameters.updateParameter *= adaptiveFactor
    if mc.parameters.updateParameter > 1.0π
        mc.parameters.updateParameter = 1.0π
    end
    return energy
end

function localSweep(mc::MonteCarlo{T}, energy::Float64) where {T<:Lattice}
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

function microcanonicalSweep!(lattice::Lattice{D,N}, rounds::Int, rng=Random.GLOBAL_RNG) where {D,N}
    basisLength = length(lattice.unitcell)
    for _ in 1:rounds
        #Deterministic first iteration
        for i in 1:basisLength
            sublatticeOrdered = range(i, length(lattice), step=basisLength)
            for site in sublatticeOrdered
                newSpinState = microcanonicalRotation(lattice, site)
                setSpin!(lattice, site, newSpinState)
            end
        end
        #Random second iteration
        for i in 1:basisLength
            sublatticeOrdered = range(i, length(lattice), step=basisLength)
            for site in sublatticeOrdered
                setSpin!(lattice, site, microcanonicalRotationRandom(lattice, site, rng))
            end
        end
    end
    return nothing
end

function microcanonicalSweep!(mc::MonteCarlo{T}) where {T<:Lattice}
    return microcanonicalSweep!(mc.lattice, mc.parameters.microcanonicalRoundsPerSweep, mc.parameters.rng)
end

"""
--------------------------------------------------------------------------------
replicaExchange functions
--------------------------------------------------------------------------------
"""

function replicaExchange!(mc::MonteCarlo{T}, energy::Float64, betas::Vector{Float64}, channelsUp, channelsDown, label) where {T<:Lattice}
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
--------------------------------------------------------------------------------
Monte Carlo Functions
--------------------------------------------------------------------------------
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

function createChannels()
    channelsUp = [fetch(@spawnat i RemoteChannel(() -> Channel{Any}(1), myid())) for i in procs()]
    channelsDown = circshift([fetch(@spawnat i RemoteChannel(() -> Channel{Any}(1), myid())) for i in procs()], -1)
    return (channelsUp, channelsDown)
end

function printStatistics!(mc::MonteCarlo{T}; replica=false) where {T<:Lattice}
    t = time()
    if mc.parameters.sweep % mc.parameters.reportInterval == 0
        #collect statistics
        totalSweeps = mc.parameters.thermalizationSweeps + mc.parameters.measurementSweeps
        progress = 100.0 * mc.parameters.sweep / totalSweeps
        thermalized = (mc.parameters.sweep >= mc.parameters.thermalizationSweeps) ? "YES" : "NO"
        sweeprate = mc.statistics.sweeps / (t - mc.statistics.initializationTime)
        sweeptime = 1.0 / sweeprate
        eta = (totalSweeps - mc.parameters.sweep) / sweeprate
        localUpdateAcceptanceRate = 100.0 * mc.statistics.acceptedLocalUpdates / mc.statistics.attemptedLocalUpdates

        #print statistics
        str = ""
        str *= @sprintf("Sweep %d / %d (%.1f%%)", mc.parameters.sweep, totalSweeps, progress)
        str *= @sprintf("\t\tETA : %s\n", Dates.format(Dates.now() + Dates.Second(round(Int64, eta)), "dd u yyyy HH:MM:SS"))
        str *= @sprintf("\t\tthermalized : %s\n", thermalized)
        str *= @sprintf("\t\tsweep rate : %.1f sweeps/s\n", sweeprate)
        str *= @sprintf("\t\tsweep duration : %.3f ms\n", sweeptime * 1000)
        str *= @sprintf("\t\tupdate acceptance rate: %.2f%%\n", localUpdateAcceptanceRate)
        str *= @sprintf("\t\tupdate parameter: %.3f \n", mc.parameters.updateParameter)
        if replica
            str *= @sprintf("\t\treplica acceptanced : %.3f \n", mc.statistics.acceptedReplicaExchanges)
            str *= @sprintf("\t\treplica attempted : %.3f \n", mc.statistics.attemptedReplicaExchanges)
        end
        str *= @sprintf("\n")
        print(str)

        #update total statistics and reset current statistics
        updateTotalStatistics!(mc.statistics)
    end
end

function sanityChecks(mc::MonteCarlo{T}, outfile::Union{String,Nothing}=nothing) where {T<:Lattice}
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
    return enableOutput
end

"""
--------------------------------------------------------------------------------
Monte Carlo run! Functions
--------------------------------------------------------------------------------
"""

function run!(mc::MonteCarlo{T}; outfile::Union{String,Nothing}=nothing) where {T<:Lattice}

    enableOutput = sanityChecks(mc, outfile)

    #init spin configuration
    initSpinConfiguration!(mc)

    #init Monte Carlo run
    totalSweeps = mc.parameters.thermalizationSweeps + mc.parameters.measurementSweeps
    energy = getEnergy(mc.lattice)

    #launch Monte Carlo run
    lastCheckpointTime = time()
    @printf("Simulation started on %s.\n\n", Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))

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
            checkpointPending = time() - lastCheckpointTime >= mc.checkpointInterval
            if checkpointPending
                writeCheckpoint!(outfile, mc)
                lastCheckpointTime = time()
                @printf("Checkpoint written on %s.\n", Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))
            end
        end
    end

    #write final checkpoint
    if enableOutput
        writeMonteCarlo!(outfile, mc)
        @printf("Checkpoint written on %s.\n", Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))
    end

    #return
    @printf("Simulation finished on %s.\n", Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))
    return nothing
end

function run!(mcs::MonteCarloAnnealing; outfile::Union{String,Nothing}=nothing)
    for (i, mc) in enumerate(mcs.MonteCarloObjects)
        if i != 1
            #If not the first simulation, copy spin configuration from the previous one
            mcs.MonteCarloObjects[i].lattice = deepcopy(mcs.MonteCarloObjects[i-1]).lattice
        end
        out = outfile === nothing ? outfile : outfile * "." * string(i - 1)
        run!(mc, outfile=out)
    end
end

function run!(mcs::MonteCarloExchange, outfile::Union{String,Nothing}=nothing)
    pmap((i) -> fetch(@spawnat i run!(mcs.MonteCarloObjects[i-1],
            mcs.betas,
            mcs.channelsUp[i-1:i], mcs.channelsDown[i-1:i],
            (outfile === nothing ? outfile : outfile * "." * string(i - 1)))),
        workers())
end

function run!(mc::MonteCarlo{T}, betas::Vector{Float64}, channelsUp::Vector{RemoteChannel{Channel{C}}},
    channelsDown::Vector{RemoteChannel{Channel{C}}},
    outfile::Union{String,Nothing}=nothing) where {T<:Lattice,C}

    enableOutput = sanityChecks(mc, outfile)
    rank = myid() - 1
    nSimulations = nworkers()
    labels = Vector{Int}(undef, floor(Int,
        (mc.parameters.thermalizationSweeps + mc.parameters.measurementSweeps) / mc.parameters.replicaExchangeRate) + 1)
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
        @printf("Simulation started on %s.\n\n", Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))
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
            energy, currentLabel = replicaExchange!(mc, energy, betas, channelsUp, channelsDown, currentLabel)
            labels[floor(Int, mc.parameters.sweep / mc.parameters.replicaExchangeRate)+1] = currentLabel
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
            printStatistics!(mc, replica=true)
        end

        #write checkpoint
        if enableOutput
            checkpointPending = time() - lastCheckpointTime >= mc.parameters.checkpointInterval
            if checkpointPending
                writeCheckpoint!(outfile, mc)
                lastCheckpointTime = time()
                rank == 1 && @printf("Checkpoint written on %s.\n", Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))
            end
        end
    end

    #write final checkpoint
    if enableOutput
        writeMonteCarlo!(outfile, mc)
        rank == 1 && @printf("Checkpoint written on %s.\n", Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))
    end

    #return
    rank == 1 && @printf("Simulation finished on %s.\n", Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))
    return nothing
end