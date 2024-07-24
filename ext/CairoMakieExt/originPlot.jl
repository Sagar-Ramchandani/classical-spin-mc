"""
--------------------------------------------------------------------------------
Building blocks for common origin plots
--------------------------------------------------------------------------------
"""

function plotArrow!(
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