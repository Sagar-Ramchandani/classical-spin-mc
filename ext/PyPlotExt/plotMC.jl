function cropXY(x, y, xLow, xHigh)
    xCrop = typeof(x)(undef, 0)
    yCrop = typeof(y)(undef, 0)
    for (i, j) in zip(x, y)
        if i > xLow && i < xHigh
            push!(xCrop, i)
            push!(yCrop, j)
        end
    end
    return xCrop, yCrop
end

function plotProperty!(ax, y; x=nothing, color=defaultColor, kwargs...)
    xAxis = x === nothing ? range(1, length(y)) : x
    yAxis = [i.val for i in y]
    Δy = [i.err for i in y]

    ax.errorbar(xAxis, yAxis, Δy, color=color, zorder=defaultZOrder, ecolor=errorColor, linewidth=0.5, elinewidth=0.5; kwargs...)
    ax.scatter(xAxis, yAxis, marker=defaultMarker, s=defaultMarkerSize, color=color, zorder=defaultZOrder + 1, edgecolors=defaultColor, linewidths=0.25)
end

function getSubplots(N::Int64; tile=false)
    if tile
        nearestSquare = floor(Int64, sqrt(N))
        additionalRows = ceil(Int64, (N - nearestSquare^2) / nearestSquare)
        return (nearestSquare + additionalRows, nearestSquare)
    else
        return (N, 1)
    end
end

function extendLimits(l, u; c=0.1)
    δ = (u - l) * c
    return (l - δ, u + δ)
end

function setTicks(ax, l::T1, u::T2) where {T1,T2}
    ax.set_yticks(round.(range(l, u, length=3), digits=2))
end

function setTicks(ax, l::Rational, u::Rational)
    nticks = 3
    ax.set_yticks(range(l, u, length=nticks), labels=[i.den == 1 ? "$(i.num)" : L"\frac{%$(i.num)}{%$(i.den)}" for i in range(l, u, length=nticks)])
end

function plotTransitionLines(axis, temperature, observables, properties)
    CvIndex = findfirst(x -> x == :specificHeat, properties)
    Cv = getindex.(observables[:, CvIndex], 1)
    transitionIndex = argmax(Cv)
    Tc = temperature[transitionIndex]
    for ax in axis
        ax.plot([Tc, Tc], [-10, 10], c=transitionColor, linestyle="-.", linewidth=0.5)
    end
end

function plotObservables!(axis, temperature::Vector{Float64}, observables::Matrix{T}, properties::Vector{Symbol};
    color=defaultColor, transition=true, markLines=true, setLimits=true, kwargs...) where {T}
    plt.minorticks_off()
    plt.xlim(0.0, maximum(temperature))
    plt.xlabel("temperature")

    transition && (:specificHeat in properties) && plotTransitionLines(axis, temperature, observables, properties)

    for (i, (ax, prop)) in enumerate(zip(axis, properties))
        y = [o.val for o in observables[:, i]]

        ax.set_ylabel(get(labelsDict, prop, prop), wrap=true)
        if setLimits
            proposedLimits = extendLimits(getLimits(y, prop)...)
            ax.set_ylim(proposedLimits)
            setTicks(ax, getTicks(y, prop)...)
        end
        markLine = getMarkerLine(y, prop)
        if markLines && !(markLine === nothing)
            if length(markLine) == 2
                ax.plot(ax.get_xlim(), markLine, c=transitionColor, zorder=defaultZOrder)
            else
                ax.plot(ax.get_xlim(), (markLine, markLine), c=transitionColor, zorder=defaultZOrder)
            end
        end
        plotProperty!(ax, observables[:, i], x=temperature, color=color; kwargs...)
    end
end

function plotBase(N::Int64)
    axisSize = getSubplots(N)
    fig, axis = subplots(axisSize..., sharex=true)
    if N == 1
        axis = [axis]
    else
        axis = permutedims(axis)
    end
    fig.align_ylabels()
    return fig, axis
end

function plotObservables(temperature::Vector{Float64}, observables::Matrix{T}, properties::Vector{Symbol};
    color=defaultColor, transition=true, markLines=true, kwargs...) where {T}
    fig, axis = plotBase(length(properties))
    plotObservables!(axis, temperature, observables, properties, color=color, transition=transition, markLines=markLines; kwargs...)
    return fig, axis
end

function plotMC!(axis, fileLocation::String, properties::Vector{Symbol}; kwargs...)
    fn = getFileNames(fileLocation)
    temperatures, β, observables = loadObservables(fn, properties)
    plotObservables!(axis, temperatures, observables, properties; kwargs...)
    return nothing
end

function plotMC(fileLocation::String, properties::Vector{Symbol}; saveLocation=nothing, kwargs...)
    fig, axis = plotBase(length(properties))
    plotMC!(axis, fileLocation::String, properties::Vector{Symbol}; kwargs...)
    !(saveLocation === nothing) && savefig("$(saveLocation)/obs.pdf", format="pdf")
    return fig, axis
end

function plotMC!(axis, fileLocations::Vector{String}, properties::Vector{Symbol}; average=true, kwargs...)
    fn = getFileNames.(fileLocations)
    data = [loadObservables(i, properties) for i in fn]
    if !average
        for (temperatures, _, observables) in data
            plotObservables!(axis, temperatures, observables, properties; kwargs...)
        end
    else
        temperatures = only(unique(getindex.(data, 1)))
        observables = getindex.(data, 3)
        meanObservables = Matrix{Measurement{Float64}}(undef, length(temperatures), length(properties))
        for i in 1:size(meanObservables, 1)
            for j in 1:size(meanObservables, 2)
                obsT = [o[i, j] for o in observables]
                meanObservables[i, j] = mean(obsT)
            end
        end
        plotObservables!(axis, temperatures, meanObservables, properties; kwargs...)
    end
end

function plotMC(fileLocations::Vector{String}, properties::Vector{Symbol}; average=true, saveLocation=nothing, kwargs...)
    fig, axis = plotBase(length(properties))
    plotMC!(axis, fileLocations::Vector{String}, properties::Vector{Symbol}; average=average, kwargs...)
    !(saveLocation === nothing) && savefig("$(saveLocation)/obs.pdf", format="pdf")
    return fig, axis
end
