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

mutable struct MonteCarloStatistics
    sweeps::Int
    attemptedLocalUpdates::Int
    acceptedLocalUpdates::Int
    attemptedReplicaExchanges::Int
    acceptedReplicaExchanges::Int
    initializationTime::Float64

    MonteCarloStatistics() = new(0, 0, 0, 0, 0, time())
end

function Base.:show(io::IO, statistics::MonteCarloStatistics)
    println(io, "MonteCarloStatistics with $(statistics.sweeps) Sweeps, initialized at $(statistics.initializationTime)")
end

mutable struct MonteCarloParameters{U<:AbstractRNG,UP<:Function}
    beta::Float64
    thermalizationSweeps::Int
    measurementSweeps::Int
    measurementRate::Int
    microcanonicalRoundsPerSweep::Int
    replicaExchangeRate::Int
    randomizeInitialConfiguration::Bool
    reportInterval::Int
    checkpointInterval::Int

    rng::U
    seed::UInt
    sweep::Int

    updateFunction::UP
    updateParameter::Float64
end

function Base.:show(io::IO, parameters::MonteCarloParameters)
    println(io, "MonteCarloParameters with β=$(parameters.beta) and update function: $(parameters.updateFunction)")
end

function MonteCarloParameters(
    beta::Float64,
    thermalizationSweeps::Int,
    measurementSweeps::Int;
    measurementRate::Int=1,
    microcanonicalRoundsPerSweep::Int=0,
    replicaExchangeRate::Int=10,
    randomizeInitialConfiguration=true,
    reportInterval::Int=round(Int, 0.05 * (thermalizationSweeps + measurementSweeps)),
    checkpointInterval::Int=3600,
    rng::U=copy(Random.GLOBAL_RNG),
    seed::UInt=rand(Random.RandomDevice(), UInt),
    sweep::Int=0,
    updateFunction::Function=sphericalUpdate,
    updateParameter::Float64=get(updateParameterDict, updateFunction, 0.0)) where {U<:AbstractRNG}

    return MonteCarloParameters(beta, thermalizationSweeps, measurementSweeps,
        measurementRate, microcanonicalRoundsPerSweep, replicaExchangeRate,
        randomizeInitialConfiguration, reportInterval, checkpointInterval,
        rng, seed, sweep, updateFunction, updateParameter)
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

function Base.:show(io::IO, mc::MonteCarlo{T}) where {T<:Lattice}
    println(io, "MonteCarlo with β=$(mc.parameters.beta) and update function: $(mc.parameters.updateFunction)")
end

function MonteCarlo(lattice::T, parameters::MonteCarloParameters,
    observables::Observables) where {T<:Lattice}
    return MonteCarlo(lattice, parameters, MonteCarloStatistics(), observables)
end

function MonteCarlo(lattice::Lattice{D,N}, parameters::Tuple{Float64,Int64,Int64}, storeAllMeasurements::Bool) where {D,N}
    pm = MonteCarloParameters(parameters...)
    obs = Observables(lattice, storeAllMeasurements)
    return MonteCarlo(lattice, pm, obs)
end

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
    initSpinConfiguration!(mc.lattice, mc.parameters.updateFunction, mc.parameters.rng)
end

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
    for i in 1:length(mc.lattice)
        site = rand(mc.parameters.rng, 1:length(mc.lattice))
        newSpinState = mc.parameters.updateFunction(mc.parameters.rng)
        energy += localUpdate(mc, site, newSpinState)
    end
    mc.statistics.sweeps += 1
    return energy
end

function localSweep(::typeof(conicalUpdate), mc::MonteCarlo{T}, energy::Float64) where {T<:Lattice}
    for i in 1:length(mc.lattice)
        site = rand(mc.parameters.rng, 1:length(mc.lattice))
        newSpinState = f(getSpin(mc.lattice, site), mc.parameters.updateParameter, mc.parameters.rng)
        energy += localUpdate(mc, site, newSpinState)
    end
    mc.statistics.sweeps += 1
    return energy
end

function localSweep(mc::MonteCarlo{T}, energy::Float64) where {T<:Lattice}
    return localSweep(mc.parameters.updateFunction, mc, energy)
end

function microcanonicalSweep!(lattice::Lattice{D,N}, rounds::Int) where {D,N}
    basisLength = length(lattice.unitcell)
    sublatticeOrdered = reduce(vcat,
        [collect(range(startIndex, length(lattice), step=basisLength)) for startIndex in 1:basisLength])
    for _ in 1:rounds
        #Deterministic first iteration
        #map((x) -> setSpin!(mc.lattice, x, microcanonicalRotation(mc.lattice, x)), sublatticeOrdered)
        for site in sublatticeOrdered
            newSpinState = microcanonicalRotation(lattice, site)
            setSpin!(mc.lattice, site, newSpinState)
        end
        #Random second iteration
        #map((x) -> setSpin!(mc.lattice, x, microcanonicalRotationRandom(mc.lattice, x, mc.parameters.rng)), sublatticeOrdered)
        for site in sublatticeOrdered
            setSpin!(mc.lattice, site, microcanonicalRotationRandom(mc.lattice, site, mc.parameters.rng))
        end
    end
    return nothing
end

function microcanonicalSweep!(mc::MonteCarlo{T}) where {T<:Lattice}
    return microcanonicalSweep!(mc.lattice, mc.parameters.microcanonicalRoundsPerSweep)
end

function printStatistics!(mc::MonteCarlo{T}, statistics::MonteCarloStatistics) where {T<:Lattice}
    t = time()
    if mc.sweep % mc.reportInterval == 0
        #collect statistics
        totalSweeps = mc.thermalizationSweeps + mc.measurementSweeps
        progress = 100.0 * mc.sweep / totalSweeps
        thermalized = (mc.sweep >= mc.thermalizationSweeps) ? "YES" : "NO"
        sweeprate = statistics.sweeps / (t - statistics.initializationTime)
        sweeptime = 1.0 / sweeprate
        eta = (totalSweeps - mc.sweep) / sweeprate
        localUpdateAcceptanceRate = 100.0 * statistics.acceptedLocalUpdates / statistics.attemptedLocalUpdates

        #print statistics
        str = ""
        str *= @sprintf("Sweep %d / %d (%.1f%%)", mc.sweep, totalSweeps, progress)
        str *= @sprintf("\t\tETA : %s\n", Dates.format(Dates.now() + Dates.Second(round(Int64, eta)), "dd u yyyy HH:MM:SS"))
        str *= @sprintf("\t\tthermalized : %s\n", thermalized)
        str *= @sprintf("\t\tsweep rate : %.1f sweeps/s\n", sweeprate)
        str *= @sprintf("\t\tsweep duration : %.3f ms\n", sweeptime * 1000)
        str *= @sprintf("\t\tupdate acceptance rate: %.2f%%\n", localUpdateAcceptanceRate)
        str *= @sprintf("\n")
        print(str)

        #reset statistics
        statistics = MonteCarloStatistics()
    end
end


function run!(mc::MonteCarlo{T}; outfile::Union{String,Nothing}=nothing) where {T<:Lattice}
    #init IO
    restrictTheta = 1.0π
    enableOutput = typeof(outfile) != Nothing
    if enableOutput
        isfile(outfile) && error("File ", outfile, " already exists. Terminating.")
    end

    if mc.microcanonicalRoundsPerSweep != 0
        for site in 1:length(mc.lattice)
            if getInteractionOnsite(mc.lattice, site) != SpinMC.InteractionMatrix(zeros(9)...)
                error("Microcanonical updates are only supported for models without on-size interactions.")
            end
        end
    end

    #init spin configuration
    if (mc.sweep == 0) && mc.randomizeInitialConfiguration
        initSpinConfiguration!(mc)
    end

    siteList = calcTriangles(mc.lattice)

    #init Monte Carlo run
    totalSweeps = mc.thermalizationSweeps + mc.measurementSweeps
    energy = getEnergy(mc.lattice)

    #launch Monte Carlo run
    lastCheckpointTime = time()
    statistics = MonteCarloStatistics()
    rank == 0 && @printf("Simulation started on %s.\n\n", Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))

    while mc.sweep < totalSweeps
        #perform local sweep
        energy = localSweep(mc, statistics, energy, restrictTheta=restrictTheta)
        #perform microcanonical sweep
        if (mc.microcanonicalRoundsPerSweep != 0) &&
           (mc.sweep % mc.microcanonicalRoundsPerSweep) == 0
            microcanonicalSweep!(mc)
        end
        #perform measurement
        if mc.sweep >= mc.thermalizationSweeps
            if mc.sweep % mc.measurementRate == 0
                performMeasurements!(mc.observables, mc.lattice, energy, siteList)
            end
        end
        #increment sweep
        statistics.sweeps += 1
        mc.sweep += 1

        #runtime statistics
        printStatistics!(mc, statistics)
        adaptiveFactor = 0.5 / (1 - (statistics.acceptedLocalUpdates / statistics.attemptedLocalUpdates))
        restrictTheta = restrictTheta * adaptiveFactor
        if restrictTheta > 1.0π
            restrictTheta = 1.0π
        end
        if (mc.sweep % mc.reportInterval) == 0
            println("adaptiveFactor is $adaptiveFactor")
            println("restrictTheta is $restrictTheta")
        end

        #write checkpoint
        if enableOutput
            checkpointPending = time() - lastCheckpointTime >= mc.checkpointInterval
            if checkpointPending
                writeMonteCarlo(outfile, mc)
                lastCheckpointTime = time()
                @printf("Checkpoint written on %s.\n", Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))
            end
        end
    end

    #write final checkpoint
    if enableOutput
        writeMonteCarlo(outfile, mc)
        @printf("Checkpoint written on %s.\n", Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))
    end

    #return
    @printf("Simulation finished on %s.\n", Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))
    return nothing
end

function replicaExchange!(mc::MonteCarlo{T}, statistics::MonteCarloStatistics, energy::Float64, allBetas::Vector{Float64}, channelsUp, channelsDown, label) where {T<:Lattice}
    #determine replica partner rank
    rank = myid() - 1
    numberWorkers = length(allBetas)
    if iseven(mc.sweep ÷ mc.replicaExchangeRate)
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
        statistics.attemptedReplicaExchanges += 1
        exchangeAccepted = false

        if iseven(rank)
            p = exp(-(allBetas[rank] - allBetas[partnerRank]) * (partnerEnergy - energy))
            exchangeAccepted = (rand(mc.rng) < min(1.0, p)) ? true : false
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

            statistics.acceptedReplicaExchanges += 1
        end
    end
    return (energy, label)
end

function anneal(mc::MonteCarlo{T}, betas::Vector{Float64}; outfile::Union{String,Nothing}=nothing) where {T<:Lattice}
    simulations = Vector{MonteCarlo}(undef, length(betas))

    for (i, beta) in enumerate(betas)
        #create one simulation for each provided beta based on the specified mc template
        simulations[i] = deepcopy(mc)
        simulations[i].beta = beta
        if i != 1
            #if this is not the first simulation, copy spin configuration from the previous one
            simulations[i].randomizeInitialConfiguration = false
            simulations[i].lattice = deepcopy(simulations[i-1].lattice)
        end
        #set outfile name for current simulation and run
        out = outfile === nothing ? outfile : outfile * "." * string(i - 1)
        run!(simulations[i], outfile=out)
    end

    return simulations
end

function run!(mc::MonteCarlo{T}, allBetas::Vector{Float64},
    channelsUp, channelsDown;
    outfile::Union{String,Expr,Nothing}=nothing, resetSpins::Bool=true) where {T<:Lattice}
    rank = myid() - 1
    nSimulations = nworkers()
    restrictTheta = 1.0π
    if rank == 1
        labels = [1]
    elseif rank == length(allBetas)
        labels = [-1]
    else
        labels = [0]
    end
    push!(mc.observables.labels, first(labels))
    println("Running replica exchanges across $(nSimulations) simulations")

    #init IO
    enableOutput = (typeof(outfile) != Nothing)
    if enableOutput
        if isfile(outfile)
            error("File $(outfile) already exists. Terminating all simulations.")
        end
    end

    if mc.microcanonicalRoundsPerSweep != 0
        for site in 1:length(mc.lattice)
            if getInteractionOnsite(mc.lattice, site) != SpinMC.InteractionMatrix(zeros(9)...)
                error("Microcanonical updates are only supported for models without on-size interactions.")
            end
        end
    end

    #init spin configuration
    if (mc.sweep == 0) && mc.randomizeInitialConfiguration
        initSpinConfiguration!(mc)
    end

    siteList = calcTriangles(mc.lattice)

    #init Monte Carlo run
    totalSweeps = mc.thermalizationSweeps + mc.measurementSweeps
    energy = getEnergy(mc.lattice)

    #launch Monte Carlo run
    lastCheckpointTime = time()
    statistics = MonteCarloStatistics()
    if rank == 1
        @printf("Simulation started on %s.\n\n", Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))
    end

    while mc.sweep < totalSweeps
        #perform local sweep
        energy = localSweep(mc, statistics, energy, restrictTheta=restrictTheta)

        #perform microcanonical sweep
        if (mc.microcanonicalRoundsPerSweep != 0) &&
           (mc.sweep % mc.microcanonicalRoundsPerSweep) == 0
            microcanonicalSweep!(mc)
        end

        #perform replica exchange
        if mc.sweep % mc.replicaExchangeRate == 0
            energy, currentLabel = replicaExchange!(mc, statistics, energy, allBetas, channelsUp, channelsDown, last(labels))
            push!(labels, currentLabel)
            push!(mc.observables.labels, currentLabel)
        end

        #perform measurement
        if mc.sweep >= mc.thermalizationSweeps
            if mc.sweep % mc.measurementRate == 0
                performMeasurements!(mc.observables, mc.lattice, energy, siteList)
            end
        end

        #increment sweep
        statistics.sweeps += 1
        mc.sweep += 1

        adaptiveFactor = 0.5 / (1 - (statistics.acceptedLocalUpdates / statistics.attemptedLocalUpdates))
        restrictTheta = restrictTheta * adaptiveFactor
        if restrictTheta > 1.0π
            restrictTheta = 1.0π
        end

        #runtime statistics
        t = time()
        if mc.sweep % mc.reportInterval == 0
            #collect statistics
            progress = 100.0 * mc.sweep / totalSweeps
            thermalized = (mc.sweep >= mc.thermalizationSweeps) ? "YES" : "NO"
            sweeprate = statistics.sweeps / (t - statistics.initializationTime)
            sweeptime = 1.0 / sweeprate
            eta = (totalSweeps - mc.sweep) / sweeprate

            localUpdateAcceptanceRate = 100.0 * statistics.acceptedLocalUpdates / statistics.attemptedLocalUpdates
            #replica statistics
            replicaExchangeAcceptanceRate = 100.0 * statistics.acceptedReplicaExchanges / statistics.attemptedReplicaExchanges
            push!(mc.observables.replicaAcceptance, replicaExchangeAcceptanceRate)
            #allLocalAppectanceRate = zeros(commSize)
            #allLocalAppectanceRate[rank + 1] = localUpdateAcceptanceRate
            #MPI.Allgather!(UBuffer(allLocalAppectanceRate, 1), MPI.COMM_WORLD)
            #allReplicaExchangeAcceptanceRate = zeros(commSize)
            #allReplicaExchangeAcceptanceRate[rank + 1] = replicaExchangeAcceptanceRate
            #MPI.Allgather!(UBuffer(allReplicaExchangeAcceptanceRate, 1), MPI.COMM_WORLD)

            #print statistics
            if rank == 1
                str = ""
                str *= @sprintf("Sweep %d / %d (%.1f%%)", mc.sweep, totalSweeps, progress)
                str *= @sprintf("\t\tETA : %s\n", Dates.format(Dates.now() + Dates.Second(round(Int64, eta)), "dd u yyyy HH:MM:SS"))
                str *= @sprintf("\t\tthermalized : %s\n", thermalized)
                str *= @sprintf("\t\tsweep rate : %.1f sweeps/s\n", sweeprate)
                str *= @sprintf("\t\tsweep duration : %.3f ms\n", sweeptime * 1000)
                str *= @sprintf("\t\treplica acceptanced : %.3f \n", statistics.acceptedReplicaExchanges)
                str *= @sprintf("\t\treplica attempted : %.3f \n", statistics.attemptedReplicaExchanges)

                #if enableMPI
                #    for n in 1:commSize
                #        str *= @sprintf("\t\tsimulation %d update acceptance rate: %.2f%%\n", n - 1, allLocalAppectanceRate[n])
                #        str *= @sprintf("\t\tsimulation %d replica exchange acceptance rate : %.2f%%\n", n - 1, allReplicaExchangeAcceptanceRate[n])
                #    end
                #else
                #    str *= @sprintf("\t\tupdate acceptance rate: %.2f%%\n", localUpdateAcceptanceRate)
                #end
                str *= @sprintf("\t\tupdate acceptance rate: %.2f%%\n", localUpdateAcceptanceRate)
                str *= @sprintf("\n")
                print(str)
            end

            #reset statistics
            statistics = MonteCarloStatistics()
        end

        #write checkpoint
        if enableOutput
            checkpointPending = time() - lastCheckpointTime >= mc.checkpointInterval
            #    enableMPI && (checkpointPending = MPIBcastBool(checkpointPending, 0, MPI.COMM_WORLD))
            if checkpointPending
                writeMonteCarlo(outfile, mc)
                lastCheckpointTime = time()
                rank == 1 && @printf("Checkpoint written on %s.\n", Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))
            end
        end
    end

    #write final checkpoint
    if enableOutput
        writeMonteCarlo(outfile, mc)
        rank == 1 && @printf("Checkpoint written on %s.\n", Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))
    end

    #return
    rank == 1 && @printf("Simulation finished on %s.\n", Dates.format(Dates.now(), "dd u yyyy HH:MM:SS"))
    return (labels)
end
