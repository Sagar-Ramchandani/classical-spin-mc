module CairoMakieExt

using CairoMakie, HDF5, StaticArrays, Statistics, StatsBase
#import ClassicalSpinMC: loadObservables, getSpins, groupSpins, copBase,
#                        plotSpins, plotSpins!, originPlot, originPlot!, gsPlot, gsPlot!,
#                        gsMultiple, originPlotMultiple,
#                        plotObservables, plotMC, plotHist, plotHistMultiple
import ClassicalSpinMC: defaultColor, defaultZOrder, sphereColor
import ClassicalSpinMC: scatterVertices3D!, plotArrow!, saveFigure!, copBase
import ClassicalSpinMC: histObservable

#include("constants.jl")
#include("io.jl")
#export loadObservables, getSpins, groupSpins

include("originPlot.jl")
export scatterVertices3D!, plotArrow!, saveFigure!, copBase
#export copBase, plotSpins, plotSpins!, originPlot, originPlot!, gsPlot, gsPlot!, gsMultiple
#include("plotMC.jl")
#export plotObservables, plotMC
include("hists.jl")
export histObservable
#export plotHist, plotHistMultiple

end
