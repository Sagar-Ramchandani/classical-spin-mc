module ClassicalSpinMC

include("Unitcell.jl")
export UnitCell, dimension, addInteraction!, setInteractionOnsite!, setField!,
       addBasisSite!, addAnisotropy!,
       resetBasis!

include("Lattice.jl")
export Lattice, size, length, getSpin, setSpin!,
       getSitePosition, getInteractionSites, getInteractionMatrices,
       getInteractionOnsite, getInteractionField

include("Observables.jl")
#Export constants
export BinnedObs, FullObs, BinnedVectorObs, FullVectorObs, VectorObs, BinnedMatrixObs,
       FullMatrixObs, MatrixObs

#Export Observables
export getSpecificHeat, AbstractObservables, Observables, performMeasurements!,
       performPostMeasurements!

include("Spin.jl")
export marsagliaSphereUpdate, sphericalUpdate, conicalUpdate,
       getEnergy, getEnergyDifference, getMagnetization, getMagnetizationPerSite,
       getCorrelation,
       calcTriangles

include("MonteCarlo.jl")
export MonteCarlo, MonteCarloStatistics, MonteCarloParameters, MonteCarloAnnealing,
       MonteCarloExchange,
       run!, anneal!, initSpinConfiguration!, localSweep, microcanonicalSweep!,
       replicaExchange!,
       localUpdate, microcanonicalRotation, microcanonicalRotationRandom, printStatistics!,
       sanityChecks,
       createChannels

include("IO.jl")
export writeUnitcell!, readUnitcell,
       writeLattice!, readLattice,
       writeMonteCarloParameters!, readMonteCarloParameters,
       writeMonteCarloStatistics!, readMonteCarloStatistics,
       writeObservables!, readObservables,
       writeMonteCarlo!, readMonteCarlo,
       save!, load #Only export for debugging

include("Plotting.jl")
export loadObservables, processObservables!, loadProcessedObservables, getSpins, groupSpins,
       copBase,
       plotSpins, plotSpins!, originPlot, originPlot!, gsPlot, gsPlot!, gsMultiple,
       originPlotMultiple,
       plotBase, plotObservables!, plotObservables, plotMC!, plotMC, plotHist,
       plotHistMultiple

using Reexport
@reexport using BinningAnalysis
@reexport using StaticArrays
end
