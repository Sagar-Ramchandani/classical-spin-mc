"""
--------------------------------------------------------------------------------
Plots building blocks
--------------------------------------------------------------------------------
"""

function plotArrow(ax, startPoint, endPoint; lengthRatio=0.05, color=defaultColor, kwargs...)
    δ = endPoint - startPoint
    ax.quiver(startPoint..., δ..., arrow_length_ratio=lengthRatio, color=color; kwargs...)
end

function scatterVertices3D!(ax, points; kwargs...)
    ax.scatter3D(getindex.(points, 1), getindex.(points, 2), getindex.(points, 3); kwargs...)
end

function plotSphere(; resolution=100, kwargs...)
    fig = plt.figure()
    axis = Axes3D(fig)
    plotSphere!(axis, resolution; kwargs...)
    return fig, axis
end

function plotSphere!(axis, resolution; color=sphereColor, kwargs...)
    N = resolution

    u = range(0, stop=2 * π, length=N)
    v = range(0, stop=π, length=N)

    x = cos.(u) * sin.(v)'
    y = sin.(u) * sin.(v)'
    z = ones(N) * cos.(v)'
    axis.plot_surface(x, y, z, alpha=0.1, shade=false, color=color; kwargs...)

    return nothing
end

function copBase(; limL=0.0, arrowL=0.0, init=(55, 45))
    fig = plt.figure()
    axis = Axes3D(fig)
    plotSphere!(axis, 100)
    axis.set_axis_off()

    (arrowL > 0.0) && map((x) -> plotArrow(axis, -x[1], x[1], color=x[2]), zip(arrowL .* XYZAxis, XYZColors))

    limFuncs = (axis.set_xlim3d, axis.set_ylim3d, axis.set_zlim3d)
    (limL > 0.0) && map((x) -> x(-limL, limL), limFuncs)

    axis.view_init(elev=init[1], azim=init[2])
    return fig, axis
end

function plotSpins!(axis, spins; arrow=false, kwargs...)
    scatterVertices3D!(axis, spins; kwargs...)
    arrow && plotArrow(axis, zeros(3), spins)
end

function plotSpins(spins; arrow=false, kwargs...)
    fig, axis = copBase()
    plotSpins!(axis, spins, arrow=arrow; kwargs...)
    return fig, axis
end

"""
--------------------------------------------------------------------------------
Functions for common origin plots
--------------------------------------------------------------------------------
"""

function originPlot!(axis, fn)
    spins = getSpins(fn)
    nSites = getNSites(fn)
    groupedSpins = groupSpins(nSites, spins)
    map((x) -> plotSpins!(axis, x[1], color=x[2]), zip(groupedSpins, XYZColors))
end

function originPlot(fn)
    fig, axis = copBase(arrowL=1.1)
    originPlot!(axis, fn)
    return fig, axis
end

function gsPlot!(axis, fileLocation; saveLocation=nothing)
    filenames = getFileNames(fileLocation)
    β = [h5read("$(fn)", "mc/parameters/beta") for fn in filenames]
    perm = sortperm(β)
    gs = last(filenames[perm])
    originPlot!(axis, gs)
    !(saveLocation === nothing) && savefig(fileLocation * "/cop.pdf", format="pdf")
    return nothing
end

function gsPlot(fileLocation; saveLocation=nothing)
    fig, axis = copBase(arrowL=1.1)
    gsPlot!(axis, fileLocation, saveLocation=saveLocation)
    return fig, axis
end

function gsMultiple(location; saveLocation=nothing)
    fn = getFolderNames(location)
    fig, axis = copBase(arrowL=1.1)
    for i in eachindex(fn)
        gsPlot!(axis, "$(fn[i])/")
    end
    !(saveLocation === nothing) && savefig(location * "/cop.pdf", format="pdf")
    return fig, axis
end