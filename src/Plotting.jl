"""
--------------------------------------------------------------------------------
Constants for all plots
--------------------------------------------------------------------------------
"""
const defaultColor = "#000000"
const defaultZOrder = 1

"""
--------------------------------------------------------------------------------
Constant for common origin plots
--------------------------------------------------------------------------------
"""
const sphereColor = "#2ec4b6"
const XYZAxis = Vector{Float64}[[1, 0, 0], [0, 1, 0], [0, 0, 1]]
const XYZColors = ["#4477AA", "#EE6677", "#228833"]
const origin = zeros(3)

"""
--------------------------------------------------------------------------------
Blank functions for extensions
--------------------------------------------------------------------------------
"""

function scatterVertices3D! end
function plotArrow! end
function saveFigure! end
function copBase end

"""
--------------------------------------------------------------------------------
Higher-level function for common origin plots
--------------------------------------------------------------------------------
"""

function plotSpins!(axis, spins; color = defaultColor, arrow = false, kwargs...)
    scatterVertices3D!(axis, spins; color = color, kwargs...)
    if arrow
        for spin in spins
            plotArrow!(axis, origin, spin, color = color)
        end
    end
end

function plotSpins(spins; arrow = false, kwargs...)
    fig, axis = copBase()
    plotSpins!(axis, spins, arrow = arrow; kwargs...)
    return fig, axis
end

function originPlot!(axis, fn; colors = XYZColors, kwargs...)
    spins = getSpins(fn)
    nSites = getNSites(fn)
    groupedSpins = groupSpins(nSites, spins)
    for (spinGroup, spinColor) in zip(groupedSpins, colors)
        plotSpins!(axis, spinGroup, color = spinColor; kwargs...)
    end
end

function originPlot(fn::String; kwargs...)
    fig, axis = copBase()
    originPlot!(axis, fn; kwargs...)
    return fig, axis
end

function originPlot(fn::Vector{String}; kwargs...)
    fig, axis = copBase()
    for ff in fn
        originPlot!(axis, ff; kwargs...)
    end
    return fig, axis
end

function gsPlot!(axis, fileLocation; saveLocation = nothing, kwargs...)
    filenames = getFileNames(fileLocation)
    β = [h5read(fn, "mc/parameters/beta") for fn in filenames]
    perm = sortperm(β)
    gs = last(filenames[perm])
    originPlot!(axis, gs; kwargs...)
    !(saveLocation === nothing) && saveFigure("$(saveLocation)/cop.pdf", axis)
    return nothing
end

function gsPlot(fileLocation::String; saveLocation = nothing, kwargs...)
    fig, axis = copBase()
    gsPlot!(axis, fileLocation, saveLocation = saveLocation; kwargs...)
    return fig, axis
end

function gsPlot(fn::Vector{String}; saveLocation = nothing, kwargs...)
    fig, axis = copBase()
    for ff in fn
        gsPlot!(axis, ff; kwargs...)
    end
    !(saveLocation === nothing) && saveFigure!(saveLocation, axis)
    return fig, axis
end

"""
--------------------------------------------------------------------------------
Higher-level functions for Monte Carlo plots.
--------------------------------------------------------------------------------
"""

function plotBase end
function plotObservables! end
function plotObservables end
function plotMC! end
function plotMC end
function plotHist end
function plotHistMultiple end
