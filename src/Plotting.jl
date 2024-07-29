"""
--------------------------------------------------------------------------------
Constants for all plots
--------------------------------------------------------------------------------
"""
const defaultColor = "#000000"
const defaultZOrder = 1

"""
--------------------------------------------------------------------------------
Constants for common origin plots
--------------------------------------------------------------------------------
"""
const sphereColor = "#2ec4b6"
const XYZAxis = Vector{Float64}[[1, 0, 0], [0, 1, 0], [0, 0, 1]]
const XYZColors = ["#4477AA", "#EE6677", "#228833"]
const origin = zeros(3)

"""
--------------------------------------------------------------------------------
Constants for MC plots
--------------------------------------------------------------------------------
"""
const defaultMarker = "o"
const defaultMarkerSize = 3
const errorColor = defaultColor
const transitionColor = "#AAAAEE"

const labelsDict = Dict{Symbol, String}(
    :energy => raw"energy $E$",
    :specificHeat => raw"specific heat $C_v$",
    :magnetization => raw"magnetization $m$"
)

const limitsDict = Dict{Symbol, Function}(
    :energy => (x) -> (floor(Int64, minimum(x)), 0.0),
    :specificHeat => (x) -> (0.0, ceil(Int64, maximum(x))),
    :magnetization => (x) -> (0.0, 2 / 3)
)

const ticksDict = Dict{Symbol, Function}(
    :energy => x -> (floor(Int64, minimum(x)), 0),
    :specificHeat => x -> (0, ceil(Int64, maximum(x))),
    :magnetization => x -> (0, 1)
)

const markerDict = Dict{Symbol, Function}(
    :energy => x -> minimum(x),
    :specificHeat => x -> 1,
    :magnetization => x -> maximum(x)
)

function getLabel(prop::Symbol)
    get(labelsDict, prop, String(prop))
end

function getLimits(observable::Vector{T}, prop::Symbol) where {T}
    defaultFunction = extrema
    return get(limitsDict, prop, defaultFunction)(observable)
end

function getTicks(observable::Vector{T}, prop::Symbol) where {T}
    defaultFunction = (x) -> round.(extrema(x), sigdigits = 3)
    return get(ticksDict, prop, defaultFunction)(observable)
end

function getMarkerLine(observable::Vector{T}, prop::Symbol) where {T}
    defaultFunction = x -> nothing
    return get(markerDict, prop, defaultFunction)(observable)
end

"""
--------------------------------------------------------------------------------
Blank functions for extensions
--------------------------------------------------------------------------------
"""

function scatterVertices3D! end
function plotArrow! end
function saveFigure! end
function copBase end

function plotBase end
function plotLine! end
function plotProperty! end
function getXLimit end
function getYLimit end
function setXLimit end
function setYLimit end
function setYLabel end
function setTicks end

function histObservable end

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

function extendLimits(l, u; c = 0.1)
    δ = (u - l) * c
    return (l - δ, u + δ)
end

function plotTransitionLines(axis, temperature, observables, properties)
    CvIndex = findfirst(x -> x == :specificHeat, properties)
    Cv = getindex.(observables[:, CvIndex], 1)
    transitionIndex = argmax(Cv)
    Tc = temperature[transitionIndex]
    for ax in axis
        plotLine!(ax, [Tc, Tc], [-10, 10], color = transitionColor, linewidth = 0.5)
        #ax.plot([Tc, Tc], [-10, 10], c = transitionColor, linestyle = "-.", linewidth = 0.5)
    end
end

function plotObservables!(axis, temperature::Vector{Float64},
        observables::Matrix{T}, properties::Vector{Symbol};
        color = defaultColor, transition = true, markLines = true, setLimits = true, kwargs...) where {T}
    maximumTemperature = maximum(temperature)

    transition && (:specificHeat in properties) &&
        plotTransitionLines(axis, temperature, observables, properties)

    for (i, (ax, prop)) in enumerate(zip(axis, properties))
        setXLimit(ax, 0.0, maximumTemperature)
        setYLabel(ax, getLabel(prop))

        y = [o.val for o in observables[:, i]]

        if setLimits
            proposedLimits = extendLimits(getLimits(y, prop)...)
            setYLimit(ax, proposedLimits...)
            setTicks(ax, getTicks(y, prop)...)
        end

        markLine = getMarkerLine(y, prop)
        if markLines && !(markLine === nothing)
            if length(markLine) == 2
                plotLine!(ax, getXLimit(ax), markLine,
                    color = transitionColor, zorder = defaultZOrder)
            else
                plotLine!(ax, getXLimit(ax), (markLine, markLine),
                    color = transitionColor, zorder = defaultZOrder)
            end
        end
        plotProperty!(ax, observables[:, i], x = temperature, color = color; kwargs...)
    end
end

function plotObservables(
        temperature::Vector{Float64}, observables::Matrix{T}, properties::Vector{Symbol};
        color = defaultColor, transition = true, markLines = true, kwargs...) where {T}
    fig, axis = plotBase(length(properties))
    plotObservables!(axis, temperature, observables, properties, color = color,
        transition = transition, markLines = markLines; kwargs...)
    return fig, axis
end

function plotMC!(axis, fileLocation::String, properties::Vector{Symbol}; kwargs...)
    fn = getFileNames(fileLocation)
    temperatures, observables = readObservable(fn, properties)
    plotObservables!(axis, temperatures, observables, properties; kwargs...)
    return nothing
end

function plotMC!(axis, fileLocations::Vector{String},
        properties::Vector{Symbol}; average = true, kwargs...)
    fn = getFileNames.(fileLocations)
    data = [readObservable(i, properties) for i in fn]
    if !average
        for (temperatures, observables) in data
            plotObservables!(axis, temperatures, observables, properties; kwargs...)
        end
    else
        temperatures = only(unique(getindex.(data, 1)))
        observables = getindex.(data, 2)
        meanObservables = Matrix{Measurement{Float64}}(
            undef, length(temperatures), length(properties))
        for i in 1:size(meanObservables, 1)
            for j in 1:size(meanObservables, 2)
                obsT = [o[i, j] for o in observables]
                meanObservables[i, j] = mean(obsT)
            end
        end
        plotObservables!(axis, temperatures, meanObservables, properties; kwargs...)
    end
end

function plotMC(
        fileLocation::String, properties::Vector{Symbol}; saveLocation = nothing, kwargs...)
    fig, axis = plotBase(length(properties))
    plotMC!(axis, fileLocation, properties; kwargs...)
    !(saveLocation === nothing) && saveFigure!(saveLocation, axis)
    return fig, axis
end

function plotMC(fileLocations::Vector{String}, properties::Vector{Symbol};
        average = true, saveLocation = nothing, kwargs...)
    fig, axis = plotBase(length(properties))
    plotMC!(axis, fileLocations, properties, average = average; kwargs...)
    !(saveLocation === nothing) && saveFigure!(saveLocation, axis)
    return fig, axis
end

"""
--------------------------------------------------------------------------------
Higher-level functions for Histograms
--------------------------------------------------------------------------------
"""

function histObservable(::Val{I}, obs::T, bins) where {I, T}
    error("A histogram method for data of dimensionality $(I) is not implemented")
end

function histObservable(obs::T, bins; kwargs...) where {T <: FullBinner}
    dimensionality = length(first(obs.x))
    histObservable(Val(dimensionality), obs.x, bins; kwargs...)
end

function histObservable(obs::Vector{T}, bins; kwargs...) where {T <: FullBinner}
    dimensionality = only(unique(length.(first.(getfield.(obs, :x)))))
    histObservable(Val(dimensionality), reduce(vcat, getfield.(obs, :x)), bins; kwargs...)
end

function histObservable(
        fn::String, obs::Symbol; bins = 100, saveLocation = nothing, kwargs...)
    h5open(fn) do f
        o = load(f["mc/observables"], String(obs))
        fig, axis = histObservable(o, bins; kwargs...)
        !(saveLocation === nothing) && saveFigure!(saveLocation, axis)
        return fig, axis
    end
end

function histObservable(
        fn::Vector{String}, obs::Symbol; bins = 100, saveLocation = nothing, kwargs...)
    o = [load(h5open(f)["mc/observables"], String(obs)) for f in fn]
    fig, axis = histObservable(o, bins; kwargs...)
    !(saveLocation === nothing) && saveFigure!(saveLocation, axis)
    return fig, axis
end