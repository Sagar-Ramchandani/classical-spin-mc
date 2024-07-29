function histBase()
    fig, axis = plt.subplots(1)
    fig.set_figwidth(600)
    fig.set_figheight(600)
    axis.set_aspect("equal")
    axis.set_axis_off()
    plt.tight_layout(pad = 0.0)

    return fig, axis
end

function histObservable(::Val{2}, obs::Vector{Vector{Float64}}, bins;
        plotrange::NTuple{2, Vector{Float64}} = (
            [-1.0, 1.0], [-1.0, 1.0]), kwargs...)
    fig, axis = histBase()
    axis.hist2d(getindex.(obs, 1), getindex.(obs, 2),
        bins = bins, range = plotrange; kwargs...)
    return fig, axis
    axis.set_xlim(plotrange[1]...)
    axis.set_ylim(plotrange[2]...)
end