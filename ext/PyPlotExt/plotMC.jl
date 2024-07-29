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

function plotLine!(axis, x, y; color = defaultColor,
        linewidth = 1.0, kwargs...)
    axis.plot(x, y, c = color, linestyle = "-.", linewidth = linewidth)
end

function plotProperty!(ax, y; x = nothing, color = defaultColor, kwargs...)
    xAxis = x === nothing ? range(1, length(y)) : x
    yAxis = [i.val for i in y]
    Δy = [i.err for i in y]

    ax.errorbar(xAxis, yAxis, Δy, color = color, zorder = defaultZOrder,
        ecolor = errorColor, linewidth = 0.5, elinewidth = 0.5; kwargs...)
    ax.scatter(xAxis, yAxis, marker = defaultMarker, s = defaultMarkerSize, color = color,
        zorder = defaultZOrder + 1, edgecolors = defaultColor, linewidths = 0.25)
end

function getSubplots(N::Int64; tile = false)
    if tile
        nearestSquare = floor(Int64, sqrt(N))
        additionalRows = ceil(Int64, (N - nearestSquare^2) / nearestSquare)
        return (nearestSquare + additionalRows, nearestSquare)
    else
        return (N, 1)
    end
end

function plotBase(N::Int64)
    axisSize = getSubplots(N)
    fig, axis = subplots(axisSize..., sharex = true)
    if N == 1
        axis = [axis]
    else
        axis = permutedims(axis)
    end
    fig.align_ylabels()
    plt.minorticks_off()
    plt.xlabel("temperature")
    return fig, axis
end

function setTicks(ax, l::T1, u::T2) where {T1, T2}
    ax.set_yticks(round.(range(l, u, length = 3), digits = 2))
end

function setTicks(ax, l::Rational, u::Rational)
    nticks = 3
    ax.set_yticks(range(l, u, length = nticks),
        labels = [i.den == 1 ? "$(i.num)" : L"\frac{%$(i.num)}{%$(i.den)}"
                  for i in range(l, u, length = nticks)])
end

getXLimit(axis) = axis.get_xlim()
getYLimit(axis) = axis.get_ylim()
setXLimit(axis, l, u) = axis.set_xlim(l, u)
setYLimit(axis, l, u) = axis.set_ylim(l, u)

setYLabel(axis, label) = axis.set_ylabel(label, wrap = true)