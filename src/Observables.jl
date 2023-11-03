using BinningAnalysis

#Defining Aliases for Observable Types
const BinnedObs = LogBinner{Float64,32,BinningAnalysis.Variance{Float64}}
const FullObs = FullBinner{Float64,Vector{Float64}}

const BinnedVectorObs = LogBinner{Vector{Float64},32,BinningAnalysis.Variance{Vector{Float64}}}
const FullVectorObs = FullBinner{Vector{Float64},Vector{Vector{Float64}}}
const VectorObs = Union{BinnedVectorObs,FullVectorObs}

const BinnedMatrixObs = LogBinner{Matrix{Float64},32,BinningAnalysis.Variance{Matrix{Float64}}}
const FullMatrixObs = FullBinner{Matrix{Float64},Vector{Matrix{Float64}}}
const MatrixObs = Union{BinnedMatrixObs,FullMatrixObs}

#function fraction(vec)
#    total=0
#    up=0
#    for i in vec
#        if i!=zero(eltype(vec))
#            total+=1
#            if i==one(eltype(vec))
#                up+=1
#            end
#        end
#    end
#    return up/total
#end

mutable struct Observables
    """
    Add fraction as a proper field and 
    possibly move labels somewhere else
    """
    energy::ErrorPropagator{Float64,32}

    magnetization::BinnedObs
    squaredMagnetization::BinnedObs
    quadMagnetization::BinnedObs
    mPlanar::BinnedObs
    z6::BinnedObs

    magnetizationVector::VectorObs
    magnetizationVectorPerSite::MatrixObs

    correlation::BinnedMatrixObs

    chirality::BinnedObs

    """
    Possibly change them to FullBinners or LogBinners
    """
    labels::Vector{Int64}
    replicaAcceptance::Vector{Float64}
end

function Observables(lattice::T, storeAll::Bool) where {T<:Lattice}
    return Observables(
        ErrorPropagator(Float64),
        LogBinner(Float64),
        LogBinner(Float64),
        LogBinner(Float64),
        LogBinner(Float64),
        LogBinner(Float64),
        storeAll ? FullBinner([zeros(Float64, 3)]) : LogBinner(zeros(Float64, 3)),
        storeAll ? FullBinner([zeros(Float64, 3, length(lattice.unitcell))]) : LogBinner(zeros(Float64, 3, length(lattice.unitcell))),
        LogBinner(zeros(Float64, lattice.length, length(lattice.unitcell))),
        LogBinner(Float64),
        Int64[],
        Float64[],
    )
end

function performMeasurements!(observables::Observables, lattice::T, energy::Float64, siteList::Vector{Vector{Int64}}) where {T<:Lattice}
    #measure energy and energy^2
    push!(observables.energy, energy / length(lattice), energy * energy / (length(lattice) * length(lattice)))

    m = getMagnetization(lattice)
    mMagnitude = norm(m)
    push!(observables.magnetizationVector, m)
    push!(observables.magnetization, mMagnitude)

    #Measuring m^2
    m2 = mMagnitude^2
    push!(observables.squaredMagnetization, m2)
    #Measuring m^4
    m4 = mMagnitude^4
    push!(observables.quadMagnetization, m4)

    #Measuring planar magnetization and z6
    mPlanar, z6 = getZ6(lattice)
    push!(observables.z6, z6)
    push!(observables.mPlanar, mPlanar)

    #measure spin correlations
    push!(observables.correlation, getCorrelation(lattice))
    push!(observables.magnetizationVectorPerSite, getMagnetizationPerSite(lattice))
    push!(observables.chirality, getChirality(lattice, siteList))
end
