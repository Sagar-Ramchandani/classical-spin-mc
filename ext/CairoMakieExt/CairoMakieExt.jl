module CairoMakieExt

using CairoMakie, HDF5, StaticArrays, Statistics, StatsBase
import ClassicalSpinMC: loadObservables, getSpins, groupSpins, copBase,
                        plotSpins, plotSpins!, originPlot, originPlot!, gsPlot, gsPlot!,
                        gsMultiple, originPlotMultiple,
                        plotObservables, plotMC, plotHist, plotHistMultiple

include("constants.jl")
include("io.jl")
export loadObservables, getSpins, groupSpins

include("originPlot.jl")
export copBase, plotSpins, plotSpins!, originPlot, originPlot!, gsPlot, gsPlot!, gsMultiple
#include("plotMC.jl")
#export plotObservables, plotMC
include("hists.jl")
export plotHist, plotHistMultiple

end
