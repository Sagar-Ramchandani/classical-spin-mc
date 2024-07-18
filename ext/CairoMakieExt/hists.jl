function getmPlanar(fn::String)
    return h5read(fn, "mc/observables/mPlanar/values")
end

function getmPlanar(fn::Vector{String})
    return reduce(hcat, getmPlanar.(fn))
end

function plotHist(
        axis, mPlanar::Matrix{Float64}, bins::Int64, plotrange::NTuple{2, Vector{Float64}})
    #h = fit(Histogram, (mPlanar[1, :], mPlanar[2, :]), nbins=bins)
    xRange = range(plotrange[1]..., length = bins)
    yRange = range(plotrange[2]..., length = bins)
    h = fit(Histogram, (mPlanar[1, :], mPlanar[2, :]), (xRange, yRange))
    hm = heatmap!(axis, h.edges[1], h.edges[2], h.weights)

    return hm
end

function plotHist(fn::Union{String, Vector{String}}, bins::Int64;
        colorbar = false,
        plotrange::NTuple{2, Vector{Float64}} = ([-1.0, 1.0], [-1.0, 1.0]),
        saveLocation::Union{String, Nothing} = nothing, kwargs...)
    fig = Figure(size = (600, 600))
    axis = colorbar ? Axis(fig[1, 1][1, 1], aspect = 1) : Axis(fig[1, 1], aspect = 1)
    hidespines!(axis)
    hidedecorations!(axis)

    hm = plotHist(axis, getmPlanar(fn), bins, plotrange)
    if colorbar
        Colorbar(fig[1, 1][2, 1], hm, vertical = false)
    end
    limits!(axis, plotrange[1]..., plotrange[2]...)

    !(saveLocation === nothing) && save("$(saveLocation)/hist2d.pdf", axis.scene)

    return fig, axis
end

function plotHistMultiple(
        fn::String, bins::Int64; saveLocation::Union{String, Nothing} = nothing, kwargs...)
    folders = getFolderNames(fn)
    files = getFileNames.(folders)
    nRuns = length(first(files))
    files = [[i[j] for i in files] for j in 1:nRuns]
    figs = []
    axises = []
    for i in eachindex(files)
        fig, axis = plotHist(files[i], bins, saveLocation = nothing; kwargs...)
        !(saveLocation === nothing) &&
            save(saveLocation * "/hist2d$(i).pdf", axis.scene, pt_per_unit = 246 / 600)
        push!(figs, fig)
        push!(axises, figs)
        println(i)
    end
    return figs, axises
end
