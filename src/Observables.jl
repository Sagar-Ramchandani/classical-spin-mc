using BinningAnalysis
using Statistics

#Defining Aliases for Observable Types
"""
    const BinnedObs
Alias for a LogBinner storing a Float64.
"""
const BinnedObs = LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
"""
    const FullObs
Alias for a FullBinner storing a Float64.
"""
const FullObs = FullBinner{Float64,Vector{Float64}}

"""
    const BinnedVectorObs
Alias for a LogBinner storing a Vector{Float64}.
"""
const BinnedVectorObs = LogBinner{Vector{Float64},32,BinningAnalysis.Variance{Vector{Float64}}}
"""
    const FullVectorObs
Alias for a FullBinner storing a Vector{Float64}.
"""
const FullVectorObs = FullBinner{Vector{Float64},Vector{Vector{Float64}}}
"""
    const VectorObs
Alias for either BinnedVectorObs or FullVectorObs.
"""
const VectorObs = Union{BinnedVectorObs,FullVectorObs}

"""
    const BinnedMatrixObs
Alias for LogBinner storing a Matrix{Float64}.
"""
const BinnedMatrixObs = LogBinner{Matrix{Float64},32,BinningAnalysis.Variance{Matrix{Float64}}}
"""
    const FullMatrixObs
Alias for FullBinner storing a Matrix{Float64}.
"""
const FullMatrixObs = FullBinner{Matrix{Float64},Vector{Matrix{Float64}}}
"""
    const MatrixObs
Alias for either BinnedMatrixObs or FullMatrixObs.
"""
const MatrixObs = Union{BinnedMatrixObs,FullMatrixObs}

"""
    abstract type AbstractObservables

This is a struct to store Observables measured during the MonteCarlo run.

!!! note In order to measure observables other than the built-in Observables type
or a subset of it, you may create a subtype of AbstractObservables 
such as MyObservables and define the following functions.

1.performMeasurements!(observables::MyObservables, lattice::T, energy::Float64) where {T<:Lattice}

2.performPostMeasurements!(observables::MyObservables, lattice::T, β::Float64) where {T<:Lattice}

Depending on the complexity of your type, you may also have to define 

1.writeObservables!(fn:H5, obs::MyObservables)

2.readObservables

if the general function cannot handle those types.
"""
abstract type AbstractObservables end

import Base: ==
function ==(a::T, b::T) where {T<:AbstractObservables}
    fields = fieldnames(T)
    status = true
    for field in fields
        status = status && (getfield(a, field) == getfield(b, field))
    end
    return status
end

"""
    mutable struct Observables
Built-in Observables type that measures simple observables and 
serves as a template for custom observable structs.

!!! note "Mandatory fields"
    In order to use replica exchange (parallel tempering)
    the labels and replicaAcceptance fields are mandatory.
"""
mutable struct Observables <: AbstractObservables
    energy::ErrorPropagator{Float64,32}
    specificHeat::Tuple{Float64,Float64}
    magnetization::BinnedObs
    magnetizationVector::VectorObs
    correlation::BinnedMatrixObs
    labels::Vector{Int64}
    replicaAcceptance::Vector{Float64}
end

#mutable struct Observables <: AbstractObservables
#"""
#Add fraction as a proper field and 
#possibly move labels somewhere else
#"""
#energy::ErrorPropagator{Float64,32}

#magnetization::BinnedObs
##squaredMagnetization::BinnedObs
##quadMagnetization::BinnedObs
##mPlanar::BinnedObs
##z6::BinnedObs

#magnetizationVector::VectorObs
##magnetizationVectorPerSite::MatrixObs

##correlation::BinnedMatrixObs

##chirality::Tuple{Vector{Vector{Int64}},BinnedObs}

#"""
#Possibly change them to FullBinners or LogBinners
#"""
#labels::Vector{Int64}
#replicaAcceptance::Vector{Float64}
#end

function Base.:show(io::IO, obs::Observables)
    println(io, "Observables object with $(length(fieldnames(typeof(obs)))) observables")
end

"""
    function Observables(lattice::T) where {T<:Lattice}
Generates an Observables object for a given Lattice object.
"""
function Observables(lattice::T) where {T<:Lattice}
    return Observables(
        ErrorPropagator(Float64),
        (0.0, 0.0),
        LogBinner(Float64),
        LogBinner(zeros(Float64, 3)),
        LogBinner(zeros(Float64, lattice.length, length(lattice.unitcell))),
        Int64[],
        Float64[],
    )
end

#function Observables(lattice::T, storeAll::Bool) where {T<:Lattice}
#return Observables(
#ErrorPropagator(Float64),
#LogBinner(Float64),
##LogBinner(Float64),
##LogBinner(Float64),
##LogBinner(Float64),
##LogBinner(Float64),
#storeAll ? FullBinner([zeros(Float64, 3)]) : LogBinner(zeros(Float64, 3)),
##storeAll ? FullBinner([zeros(Float64, 3, length(lattice.unitcell))]) : LogBinner(zeros(Float64, 3, length(lattice.unitcell))),
##LogBinner(zeros(Float64, lattice.length, length(lattice.unitcell))),
##(calcTriangles(lattice), LogBinner(Float64)),
#Int64[],
#Float64[],
#)
#end

#Review this
function getSpecificHeat(obs::Observables, β::Float64, N::Int64)
    k = β * β * N
    c(e) = k * (e[2] - e[1] * e[1])
    ∇c(e) = [-2.0 * e[1] * k, k]
    heat = mean(obs.energy, c)
    reliableEnergyLevel = BinningAnalysis._reliable_level(obs.energy)
    dheat = sqrt(abs(var(obs.energy, ∇c, reliableEnergyLevel)) / obs.energy.count[reliableEnergyLevel])
    return (heat, dheat)
end

function getFraction(vec)
    total = 0
    up = 0
    for i in vec
        if i != zero(eltype(vec))
            total += 1
            if i == one(eltype(vec))
                up += 1
            end
        end
    end
    return up / total
end

#Method for generic computation of observables
#function performMeasurements!(obs::O, lattice::T, energy::Float64) where {O<:AbstractObservables,T<:Lattice}
#fields = fieldnames(O)
#for field in fields
#push!(getfield(obs, field), getfield(Main, Symbol("get" * String(field)))(lattice))
#end
#return nothing
#end

"""
    function performMeasurements!(observables::Observables, lattice::T, energy::Float64) where {T<:Lattice}
Performs measurements and is called during measurement sweeps. 

!!! note "Energy measurement"
    Since the current energy of the spin configuration is calculated as part of the 
    Monte Carlo sweep, it is passed on directly.
"""
function performMeasurements!(observables::Observables, lattice::T, energy::Float64) where {T<:Lattice}
    #measure energy and energy^2
    push!(observables.energy, energy / length(lattice), energy * energy / (length(lattice) * length(lattice)))

    m = getMagnetization(lattice)
    mMagnitude = norm(m)
    push!(observables.magnetizationVector, m)
    push!(observables.magnetization, mMagnitude)

    #measuring correlation
    push!(observables.correlation, getCorrelation(lattice))
end

##Specialized method for built-in observables type
#function performMeasurements!(observables::Observables, lattice::T, energy::Float64) where {T<:Lattice}
##measure energy and energy^2
##siteList = first(observables.chirality)
#push!(observables.energy, energy / length(lattice), energy * energy / (length(lattice) * length(lattice)))

#m = getMagnetization(lattice)
#mMagnitude = norm(m)
#push!(observables.magnetizationVector, m)
#push!(observables.magnetization, mMagnitude)

##Measuring m^2
##m2 = mMagnitude^2
##push!(observables.squaredMagnetization, m2)
##Measuring m^4
##m4 = mMagnitude^4
##push!(observables.quadMagnetization, m4)

##Measuring planar magnetization and z6
##mPlanar, z6 = getZ6(lattice)
##push!(observables.z6, z6)
##push!(observables.mPlanar, mPlanar)

##measure spin correlations
##push!(observables.correlation, getCorrelation(lattice))
##push!(observables.magnetizationVectorPerSite, getMagnetizationPerSite(lattice))
##push!(last(observables.chirality), getChirality(lattice, siteList))
#end

"""
    function performPostMeasurements!(observables::Observables, lattice::T, β::Float64) where {T<:Lattice}
Performs measurements and is called after all Monte Carlo sweeps are performed. 
It is intended to do any post measurements based on the collect observables or final state of the system. 
"""
function performPostMeasurements!(observables::Observables, lattice::T, β::Float64) where {T<:Lattice}
    observables.specificHeat = getSpecificHeat(observables, β, size(lattice))
end