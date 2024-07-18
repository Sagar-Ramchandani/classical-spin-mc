"""
--------------------------------------------------------------------------------
Plots building blocks
--------------------------------------------------------------------------------
"""

function plotArrow(
        ax, startPoint, endPoint; lengthRatio = 1, color = defaultColor, kwargs...)
    δ = endPoint - startPoint
    arrows!(ax, [[i] for i in startPoint]..., [[i] for i in δ]...,
        color = color, arrowsize = [0.1, 0.1, 0.3])
end

function scatterVertices3D!(
        ax, points; color = defaultColor, markersize = 5, distanceShade = false, kwargs...)
    if length(points) > 200
        points = rand(points, 200)
    end
    if distanceShade
        camerapos = camera(ax.scene).eyeposition.val
        distances = sqrt.([sum((i - camerapos) .* (i - camerapos)) for i in points])
        minD, maxD = extrema(distances)
        distances = (distances .- minD) / (maxD - minD)
        distances ./= 2
        color = [RGBAf(color.r, color.g, color.b, 1 - i) for i in distances]
    end
    scatter!(ax, getindex.(points, 1), getindex.(points, 2), getindex.(points, 3),
        markersize = markersize, color = color; kwargs...)
end

function plotSphere(; resolution = 100, kwargs...)
    fig = Figure(backgroundcolor = :transparent)
    axis = Axis3(fig)
    plotSphere!(axis, resolution; kwargs...)
    return fig, axis
end

function plotSphere!(axis, resolution; color = [sphereColor], kwargs...)
    N = resolution

    u = range(0, stop = 2 * π, length = N)
    v = range(0, stop = π, length = N)

    x = cos.(u) * sin.(v)'
    y = sin.(u) * sin.(v)'
    z = ones(N) * cos.(v)'
    surface!(axis, x, y, z, colormap = color, shading = NoShading,
        alpha = 0.1, rasterize = true; kwargs...)

    return nothing
end

function copBase(; arrowL = 0.0, limL = 1, init = (acos(1 / sqrt(3)), π / 4))
    #fig = Figure(backgroundcolor=:transparent)
    fig = Figure(size = (600, 600))
    axis = Axis3(fig[1, 1], aspect = :equal, elevation = init[1],
        azimuth = init[2], protrusions = 0, viewmode = :fit)
    hidedecorations!(axis)
    hidespines!(axis)
    plotSphere!(axis, 100)

    (arrowL > 0.0) && map((x) -> plotArrow(axis, -x[1], x[1], color = x[2]),
        zip(arrowL .* XYZAxis, XYZColors))

    (limL > 0.0) && limits!(axis, -limL, limL, -limL, limL, -limL, limL)

    return fig, axis
end

function plotSpins!(axis, spins; color = defaultColor, arrow = false, kwargs...)
    scatterVertices3D!(axis, spins; color = color, kwargs...)
    if arrow
        for s in spins
            plotArrow(axis, zeros(3), s, color = color)
        end
    end
end

function plotSpins(spins; arrow = false, kwargs...)
    fig, axis = copBase()
    plotSpins!(axis, spins, arrow = arrow; kwargs...)
    return fig, axis
end

"""
--------------------------------------------------------------------------------
Functions for common origin plots
--------------------------------------------------------------------------------
"""

function originPlot!(axis, fn; kwargs...)
    spins = getSpins(fn)
    nSites = getNSites(fn)
    groupedSpins = groupSpins(nSites, spins)
    map((x) -> plotSpins!(axis, x[1], color = x[2]; kwargs...),
        zip(groupedSpins, XYZColors))
end

function originPlot(fn)
    fig, axis = copBase(arrowL = 1.2)
    originPlot!(axis, fn)
    return fig, axis
end

function gsPlot!(axis, fileLocation; saveLocation = nothing, kwargs...)
    filenames = getFileNames(fileLocation)
    β = [h5read("$(fn)", "mc/parameters/beta") for fn in filenames]
    perm = sortperm(β)
    gs = last(filenames[perm])
    originPlot!(axis, gs; kwargs...)
    !(saveLocation === nothing) && save(fileLocation * "/cop.pdf", axis.scene)
    return nothing
end

function gsPlot(fileLocation; saveLocation = nothing)
    fig, axis = copBase(arrowL = 1.5)
    gsPlot!(axis, fileLocation, saveLocation = saveLocation)
    return fig, axis
end

function gsMultiple(location; saveLocation = nothing)
    fn = getFolderNames(location)
    fig, axis = copBase(arrowL = 1.5)
    for i in eachindex(fn)
        gsPlot!(axis, "$(fn[i])/")
    end
    !(saveLocation === nothing) && save(location * "/cop.pdf", axis.scene)
    return fig, axis
end

function originPlotMultiple(location; saveLocation = nothing)
    fn = getFolderNames(location)
    filesPerFolder = getFileNames.(fn)
    nRuns = length(first(filesPerFolder))
    filesPerFolder = [[i[j] for i in filesPerFolder] for j in 1:nRuns]
    figs = []
    axises = []
    for i in eachindex(filesPerFolder)
        fig, axis = copBase(arrowL = 1.5)
        for j in filesPerFolder[i]
            originPlot!(axis, j)
        end
        push!(figs, fig)
        push!(axises, figs)
        !(saveLocation === nothing) &&
            save(saveLocation * "/cop$(i).pdf", axis.scene, pt_per_unit = 246 / 600)
        println(i)
    end
    return figs, axises
end
