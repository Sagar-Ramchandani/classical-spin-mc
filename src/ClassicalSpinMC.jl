module ClassicalSpinMC

include("Unitcell.jl")
export UnitCell, addInteraction!, setInteractionOnsite!, setField!, addBasisSite!

include("Lattice.jl")
export Lattice, size, length, getSpin, setSpin!,
    getSitePosition, getInteractionSites, getInteractionMatrices,
    getInteractionOnsite, getInteractionField

include("Observables.jl")
export Observables

include("Spin.jl")
export marsagliaSphereUpdate, sphericalUpdate, conicalUpdate,
    getEnergy, getMagnetization, getMagnetizationPerSite, getCorrelation

include("MonteCarlo.jl")
export MonteCarlo, MonteCarloStatistics, MonteCarloParameters, run!, anneal, initSpinConfiguration!, localSweep, microcanonicalSweep!, replicaExchange!

include("IO.jl")
export writeUnitcell!, readUnitcell,
    writeLattice!, readLattice,
    writeMonteCarlo!, readMonteCarlo

using Reexport
@reexport using BinningAnalysis
@reexport using StaticArrays
end
