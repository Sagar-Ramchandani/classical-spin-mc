function histBase()
    fig = Figure(size = (600, 600))
    axis = Axis(fig[1, 1], aspect = 1)
    hidespines!(axis)
    hidedecorations!(axis)

    return fig, axis
end

function histObservable(::Val{2}, obs::Vector{Vector{Float64}}, bins;
        plotrange::NTuple{2, Vector{Float64}} = (
            [-1.0, 1.0], [-1.0, 1.0]), kwargs...)
    fig, axis = histBase()
    xRange = range(plotrange[1]..., length = bins)
    yRange = range(plotrange[2]..., length = bins)
    h = fit(Histogram, (getindex.(obs, 1), getindex.(obs, 2)), (xRange, yRange))
    hm = heatmap!(axis, h.edges[1], h.edges[2], h.weights; kwargs...)

    limits!(axis, plotrange[1]..., plotrange[2]...)

    return fig, axis
end