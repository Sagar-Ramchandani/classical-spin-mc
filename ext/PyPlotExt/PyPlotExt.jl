module PyPlotExt

using PyCall, PyPlot, HDF5, StaticArrays, Statistics, Measurements

import ClassicalSpinMC: defaultColor, defaultZOrder, sphereColor, errorColor,
                        transitionColor, defaultMarker, defaultMarkerSize
import ClassicalSpinMC: scatterVertices3D!, plotArrow!, copBase, saveFigure!
import ClassicalSpinMC: histObservable
import ClassicalSpinMC: extendLimits, plotBase, setTicks, plotLine!, plotProperty!,
                        getXLimit, getYLimit, setXLimit, setYLimit, setYLabel, setTicks

function saveFigure!(location, axis; format = "pdf")
    savefig(location, format = format)
end

include("originPlot.jl")
export scatterVertices3D!, plotArrow!, copBase, saveFigure!

include("plotMC.jl")
export setTicks, plotBase, plotLine!, plotProperty!, getXLimit, getYLimit, setXLimit,
       setYLimit, setYLabel

include("hists.jl")
export histObservable

end
