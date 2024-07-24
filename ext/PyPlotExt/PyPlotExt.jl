module PyPlotExt

using PyCall, PyPlot, HDF5, StaticArrays, Statistics, Measurements
#import ClassicalSpinMC: readObservable, processObservables!, loadProcessedObservables,
#                        getSpins, groupSpins,
#                        plotSpins, plotSpins!, originPlot, originPlot!, gsPlot, gsPlot!,
#                        gsMultiple,
#                        plotBase, plotObservables!, plotObservables, plotMC!, plotMC,
#                        plotHist

import ClassicalSpinMC: defaultColor, defaultZOrder, sphereColor
import ClassicalSpinMC: scatterVertices3D!, plotArrow!, copBase, saveFigure!
import ClassicalSpinMC: histObservable

function saveFigure!(location, axis; format = "pdf")
    savefig(location, format = format)
end

include("originPlot.jl")
export scatterVertices3D!, plotArrow!, copBase, saveFigure!

#include("plotMC.jl")
#export plotBase, plotObservables!, plotObservables, plotMC!, plotMC

include("hists.jl")
export histObservable
#export plotHist
end
