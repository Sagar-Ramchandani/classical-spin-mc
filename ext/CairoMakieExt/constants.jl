"""
--------------------------------------------------------------------------------
Constants for all plots
--------------------------------------------------------------------------------
"""
const defaultColor = "#000000"
const defaultZOrder = 1
const H5 = Union{HDF5.File, HDF5.Group}

"""
------------------------------------------------------------------------------
Constants for MC plots
------------------------------------------------------------------------------
"""
const defaultMarker = "o"
const defaultMarkerSize = 10
const errorColor = defaultColor
const transitionColor = "#AAAAEE"

const labelsDict = Dict(
    :energy => raw"energy $E$",
    :specificHeat => raw"specific heat $C_v$",
    :magnetization => raw"magnetization $m$",
    :m2 => raw"$m^2 ",
    :m4 => raw"$m^4 ",
    :mPlanar => raw"planar magnetization $m_p $",
    :z6 => raw"order parameter $z_6$",
    :chirality => raw"chirality $\chi ",
    :binder => raw"binder cumulant $U_L$"
)

function getLimits(observable::Vector{T}, prop::Symbol) where {T}
    defaultFunction = extrema
    limitsDict = Dict{Symbol, Any}(
        :energy => (x) -> (minimum(x), 0.0),
        :specificHeat => (x) -> (0.0, maximum(x)),
        :magnetization => (x) -> (0.0, 2 / 3),
        :m2 => (x) -> (0.0, (1)^2),
        :m4 => (x) -> (0.0, (1)^4),
        :mPlanar => (x) -> (0.0, maximum(x)),
        #:z6 => x -> (mean(x) > 0) ? (0.0, maximum(x)) : (minimum(x), 0.0),
        :z6 => x -> (mean(x) > 0) ? (0.0, 1.0) : (-1.0, 0.0),
        :chirality => (x) -> (-2 * sqrt(6) / 9, 2 * sqrt(6) / 9),
        :fraction => (x) -> (0, 1),
        :replicaAcceptance => (x) -> (0, maximum(x))
    )
    return get(limitsDict, prop, defaultFunction)(observable)
end

function getTicks(observable::Vector{T}, prop::Symbol) where {T}
    defaultFunction = (x) -> round.(extrema(x), sigdigits = 3)
    ticksDict = Dict(
        :energy => x -> (minimum(x), 0),
        :specificHeat => x -> (0, maximum(x)),
        :magnetization => x -> (0 // 1, 2 // 3),
        :m2 => x -> (0, (1)^2),
        :m4 => x -> (0, (1)^4),
        :mPlanar => x -> (0, maximum(x)),
        #:z6 => x -> ((mean(x) > 0) ? (0.0, maximum(x)) : (minimum(x), 0.0)),
        :z6 => x -> (mean(x) > 0) ? (0.0, 1.0) : (-1.0, 0.0),
        :chirality => x -> (-2 * sqrt(6) / 9, 2 * sqrt(6) / 9),
        :fraction => x -> (0, 1),
        :replicaAcceptance => x -> (0, maximum(x))
    )
    return get(ticksDict, prop, defaultFunction)(observable)
end

function getMarkerLine(observable::Vector{T}, prop::Symbol) where {T}
    defaultFunction = x -> nothing
    markerDict = Dict(
        :specificHeat => x -> 1,
        :magnetization => x -> maximum(x),
        :mPlanar => x -> maximum(x),
        :z6 => x -> 0,
        :chirality => x -> 0,
        :energy => x -> minimum(x),
        :fraction => x -> (1, 0)
    )
    return get(markerDict, prop, defaultFunction)(observable)
end

"""
--------------------------------------------------------------------------------
Constants for common origin plots
--------------------------------------------------------------------------------
"""
const sphereColor = colorant"#2EC4B6FF"
const XYZAxis = Vector{Float64}[[1, 0, 0], [0, 1, 0], [0, 0, 1]]
const XYZColors = [colorant"#4477AAFF", colorant"#EE6677FF", colorant"#228833FF"]
