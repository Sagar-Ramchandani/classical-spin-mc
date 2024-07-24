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

function plotHist(fn::String, bins::Int64;
        colorbar = false,
        plotrange::NTuple{2, NTuple{2, Float64}} = ((-1.0, 1.0), (-1.0, 1.0)),
        saveLocation::Union{String, Nothing} = nothing,
        kwargs...)
    fig, axis = plt.subplots(1)
    mPlanar = Vector{Matrix{Float64}}(undef, 1)
    h5open(fn, "r") do f
        mPlanar = read(f["mc/observables/mPlanar/values"])
    end
    plt.hist2d(mPlanar[1, :], mPlanar[2, :], bins = bins, range = plotrange; kwargs...)

    axis.set_aspect("equal")
    axis.set_axis_off()
    plt.tight_layout(pad = 0.0)
    colorbar && plt.colorbar()

    !(saveLocation === nothing) && savefig("$(saveLocation)/hist2d.pdf", format = "pdf",
        bbox_inches = "tight", pad_inches = 0.0)

    return fig, axis
end
