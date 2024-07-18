module PyPlotExt

using PyCall, PyPlot, HDF5, StaticArrays, Statistics, Measurements
import ClassicalSpinMC: loadObservables, processObservables!, loadProcessedObservables,
                        getSpins, groupSpins,
                        plotSpins, plotSpins!, originPlot, originPlot!, gsPlot, gsPlot!,
                        gsMultiple,
                        plotBase, plotObservables!, plotObservables, plotMC!, plotMC,
                        plotHist

using3D()
include("constants.jl")
#include("config.jl")
include("io.jl")
export loadObservables, processObservables!, loadProcessedObservables, getSpins, groupSpins

include("originPlot.jl")
export plotSpins, plotSpins!, originPlot, originPlot!, gsPlot, gsPlot!, gsMultiple
include("plotMC.jl")
export plotBase, plotObservables!, plotObservables, plotMC!, plotMC
include("hists.jl")
export plotHist

end #module
