module PyPlotExt

using PyCall, PyPlot, HDF5, StaticArrays, Statistics
import ClassicalSpinMC: loadObservables, getSpins, groupSpins,
    plotSpins, plotSpins!, originPlot, originPlot!, gsPlot, gsPlot!, gsMultiple,
    plotObservables, plotMC


using3D()
include("constants.jl")
#include("config.jl")
include("io.jl")
export loadObservables, getSpins, groupSpins

include("originPlot.jl")
export plotSpins, plotSpins!, originPlot, originPlot!, gsPlot, gsPlot!, gsMultiple
include("plotMC.jl")
export plotObservables, plotMC

end #module
