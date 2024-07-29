function plotLine!(axis, x, y; color = defaultColor, linewidth = 1.0, kwargs...)
    lines!(axis, collect(x), collect(y), color = color, linewidth = linewidth)
end

function plotProperty!(ax, y; x = nothing, color = defaultColor, kwargs...)
    xAxis = x === nothing ? range(1, length(y)) : x
    yAxis = [i.val for i in y]
    Δy = [i.err for i in y]

    lines!(ax, xAxis, yAxis, color = color)
    errorbars!(ax, xAxis, yAxis, Δy, color = errorColor, linewidth = 0.5)
    scatter!(ax, xAxis, yAxis, markersize = defaultMarkerSize, color = color)
end

function plotBase(N::Int64)
    fig = Figure()
    axis = [Axis(fig[i, 1], xgridvisible = false,
                ygridvisible = false, xlabel = latexstring("temperature")) for i in 1:N]
    linkxaxes!(axis...)
    yspace = maximum(tight_yticklabel_spacing!, axis)

    for ax in axis
        ax.yticklabelspace = yspace
    end

    return fig, axis
end

function setTicks(ax, l, u)
    ax.yticks = round.(range(l, u, length = 3), digits = 2)
end

function getXLimit(axis)
    return axis.limits.val[1]
end

function getYLimit(axis)
    return axis.limits.val[2]
end

setXLimit(axis, l, u) = xlims!(axis, l, u)
setYLimit(axis, l, u) = ylims!(axis, l, u)
function setYLabel(axis, label)
    axis.ylabel = latexstring(label)
end