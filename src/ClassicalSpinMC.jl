module ClassicalSpinMC

include("Unitcell.jl")
export UnitCell, addInteraction!, setInteractionOnsite!, setField!, addBasisSite!, addAnisotropy!

include("Lattice.jl")
export Lattice, size, length, getSpin, setSpin!,
    getSitePosition, getInteractionSites, getInteractionMatrices,
    getInteractionOnsite, getInteractionField

include("Observables.jl")
export AbstractObservables, Observables, performMeasurements!

include("Spin.jl")
export marsagliaSphereUpdate, sphericalUpdate, conicalUpdate,
    getEnergy, getEnergyDifference, getMagnetization, getMagnetizationPerSite, getCorrelation,
    calcTriangles

include("MonteCarlo.jl")
export MonteCarlo, MonteCarloStatistics, MonteCarloParameters, MonteCarloAnnealing, MonteCarloExchange,
    run!, anneal, initSpinConfiguration!, localSweep, microcanonicalSweep!, replicaExchange!,
    localUpdate, microcanonicalRotation, microcanonicalRotationRandom, printStatistics!, sanityChecks,
    createChannels

include("IO.jl")
export writeUnitcell!, readUnitcell,
    writeLattice!, readLattice,
    writeMonteCarloParameters!, readMonteCarloParameters,
    writeObservables!, readObservables,
    writeMonteCarlo!, readMonteCarlo,
    save!, load #Only export for debugging

using Reexport
@reexport using BinningAnalysis
@reexport using StaticArrays
end
