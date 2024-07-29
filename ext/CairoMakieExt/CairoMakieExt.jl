module CairoMakieExt

using CairoMakie, HDF5, StaticArrays, Statistics, StatsBase, LaTeXStrings

import ClassicalSpinMC: defaultColor, defaultZOrder, sphereColor, errorColor,
                        transitionColor, defaultMarker, defaultMarkerSize
import ClassicalSpinMC: scatterVertices3D!, plotArrow!, saveFigure!, copBase
import ClassicalSpinMC: histObservable
import ClassicalSpinMC: extendLimits, plotBase, setTicks, plotLine!, plotProperty!,
                        getXLimit, getYLimit, setXLimit, setYLimit, setYLabel, setTicks

function saveFigure!(location, axis)
    save(location, axis.scene)
end

include("originPlot.jl")
export scatterVertices3D!, plotArrow!, saveFigure!, copBase

include("plotMC.jl")
export setTicks, plotBase, plotLine!, plotProperty!, getXLimit, getYLimit, setXLimit,
       setYLimit, setYLabel

include("hists.jl")
export histObservable

end
