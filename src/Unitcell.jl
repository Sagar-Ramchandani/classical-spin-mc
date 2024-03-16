using StaticArrays

"""
    mutable struct UnitCell

This is used to represent a UnitCell struct. 
For further computation it must first be converted into a Lattice struct.
"""
@kwdef mutable struct UnitCell{D}
    primitive::NTuple{D,SVector{D,Float64}}
    basis::Vector{SVector{D,Float64}} = Vector{SVector{D,Float64}}(undef, 0)
    #interactions specified as (basis1=>basis2,offsetPrimitive,M)
    interactions::Vector{Tuple{Pair{Int,Int},NTuple{D,Int},SMatrix{3,3,Float64,9}}} =
        Vector{Tuple{Pair{Int,Int},NTuple{D,Int},SMatrix{3,3,Float64,9}}}(undef, 0)
    interactionsOnsite::Vector{SMatrix{3,3,Float64,9}} = Vector{SMatrix{3,3,Float64,9}}(undef, 0)
    interactionsField::Vector{SVector{3,Float64}} = Vector{SVector{3,Float64}}(undef, 0)
    anisotropyFunction::Function = identity
    anisotropyParameters::Vector{Float64} = [zero(Float64)]
end

"""
--------------------------------------------------------------------------------
Constructor functions
--------------------------------------------------------------------------------
"""

"""
    function UnitCell(primitive::Vararg{SVector{D,Float64},D})
UnitCell constructor that accepts D number of D dimensional SVectors as input to use as primitives. 
"""
function UnitCell(primitive::Vararg{SVector{D,Float64},D}) where {D}
    return UnitCell{D}(primitive=primitive)
end

"""
    function UnitCell(primitive::Vararg{NTuple{D,Float64},D})
UnitCell constructor that accepts D number of D dimensional NTuples for primitives that get converted to SVectors.
"""
function UnitCell(primitive::Vararg{NTuple{D,Float64},D}) where {D}
    return UnitCell(SVector.(primitive)...)
end

"""
    function addAnisotropy!(unitcell::UnitCell{D}, func::F, funcParameters::Vector{Float64}) where {D,F<:Function)
Adds an anisotropy function to the UnitCell definition along with any parameters that need to be passed 
along to the anisotropy function. The default anisotropy function is identity.
"""
function addAnisotropy!(unitcell::UnitCell{D}, func::F, funcParameters::Vector{Float64}) where {D,F<:Function}
    unitcell.anisotropyFunction = func
    unitcell.anisotropyParameters = funcParameters
end

"""
    function dimension(unitcell::UnitCell{D}) where {D}
Utility function to get the dimensionality of a UnitCell struct.
"""
function dimension(unitcell::UnitCell{D}) where {D}
    return first(typeof(unitcell).parameters)
end

"""
--------------------------------------------------------------------------------
addInteraction! functions
--------------------------------------------------------------------------------
"""

"""
    function addInteraction!(unitcell::UnitCell{D}, b::Pair{Int,Int}, M::SMatrix{3,3,Float64,9}, offset::NTuple{D,Int}=Tuple(zeros(Int, D))) where {D}
Adds an interaction between `spin_1` located at basis site `b_1` of the given `unitcell` and `spin_2` 
at basis site `b_2` in a unit cell that is offset by `offset` lattice vectors. 
The exchange energy is calculated as `spin_1'M.spin_2`. 

Please note that local interactions that involve the same site will be automatically 
added as an on site interaction using `setInteractionOnSite!(...)`.
"""
function addInteraction!(unitcell::UnitCell{D}, b::Pair{Int,Int}, M::SMatrix{3,3,Float64,9}, offset::NTuple{D,Int}=Tuple(zeros(Int, D))) where {D}
    if (b.first == b.second) && (offset == Tuple(zeros(Int, D)))
        @warn "Local Interaction detected, using setInteractionOnSite!() instead"
        setInteractionOnsite!(unitcell, b.first, M)
    end
    push!(unitcell.interactions, (b, offset, M))
    return nothing
end

"""
    function addInteraction!(unitcell::UnitCell{D}, b1::Int, b2::Int, M::SMatrix{3,3,Float64,9}, offset::NTuple{D,Int}=Tuple(zeros(Int, D))) where {D}
Convenience function that creates the `b1` and `b2` pair automatically.
"""
function addInteraction!(unitcell::UnitCell{D}, b1::Int, b2::Int, M::SMatrix{3,3,Float64,9}, offset::NTuple{D,Int}=Tuple(zeros(Int, D))) where {D}
    addInteraction!(unitcell, b1 => b2, M, offset)
end

"""
    function addInteraction!(unitcell::UnitCell{D}, b1::Int, b2s::Vector{Int}, M::SMatrix{3,3,Float64,9}, offsets::Vector{NTuple{D,Int}}) where {D}
Convenience function that adds the same interaction for 
a vector of `b2` and `offset` values pairwise.
"""
function addInteraction!(unitcell::UnitCell{D}, b1::Int, b2s::Vector{Int}, M::SMatrix{3,3,Float64,9}, offsets::Vector{NTuple{D,Int}}) where {D}
    for (b2, offset) in zip(b2s, offsets)
        addInteraction!(unitcell, b1 => b2, M, offset)
    end
    return nothing
end

"""
    function addInteraction!(unitcell::UnitCell{D}, b1::Int, b2::Int, M::Matrix{Float64}, offset::NTuple{D,Int}=Tuple(zeros(Int, D))) where {D}
Convenience function that converts to SMatrix and creates the `b1` and `b2` pair automatically.
"""
function addInteraction!(unitcell::UnitCell{D}, b1::Int, b2::Int, M::Matrix{Float64}, offset::NTuple{D,Int}=Tuple(zeros(Int, D))) where {D}
    staticM = SMatrix{size(M)...}(M)
    addInteraction!(unitcell, b1 => b2, staticM, offset)
end

"""
--------------------------------------------------------------------------------
setInteractionOnSite! functions
--------------------------------------------------------------------------------
"""

"""
    function setInteractionOnsite!(unitcell::UnitCell{D}, b::Int, M::SMatrix{3,3,Float64,9}) where {D}
Sets the self-interaction for a spin located at basis site 'b' of the given 'unitcell'. 
The exchange energy is calculated as spin'.M.spin
"""
function setInteractionOnsite!(unitcell::UnitCell{D}, b::Int, M::SMatrix{3,3,Float64,9}) where {D}
    unitcell.interactionsOnsite[b] = M
    return nothing
end

"""
    function setInteractionOnsite!(unitcell::UnitCell{D}, M::Vector{SMatrix{3,3,Float64,9}}) where {D}
Sets the self-interaction for all spins of the given 'unitcell'. 

The exchange energy for each spin is calculated as spins'[i].M[i].spins[i]
"""
function setInteractionOnsite!(unitcell::UnitCell{D}, M::Vector{SMatrix{3,3,Float64,9}}) where {D}
    unitcell.interactionsOnsite = M
    return nothing
end

"""
--------------------------------------------------------------------------------
setField! functions
--------------------------------------------------------------------------------
"""

"""
    function setField!(unitcell::UnitCell{D}, b::Int, B::SVector{3,Float64}) where {D}
Sets the magnetic field for a spin at the basis site 'b'.
"""
function setField!(unitcell::UnitCell{D}, b::Int, B::SVector{3,Float64}) where {D}
    unitcell.interactionsField[b] = B
    return nothing
end

"""
    function setField!(unitcell::UnitCell{D}, B::Vector{SVector{3,Float64}}) where {D}
Sets the magnetic field for all spins of the given unitcell. 
"""
function setField!(unitcell::UnitCell{D}, B::Vector{SVector{3,Float64}}) where {D}
    unitcell.interactionsField = B
    return nothing
end

"""
--------------------------------------------------------------------------------
Basis functions
--------------------------------------------------------------------------------
"""

"""
    function addBasisSite!(unitcell::UnitCell{D}, position::SVector{D,Float64}) where {D}
Adds a new basis site to the given `unitcell` at the given position with no self-interaction and no magnetic field.
"""
function addBasisSite!(unitcell::UnitCell{D}, position::SVector{D,Float64}) where {D}
    push!(unitcell.basis, position)
    push!(unitcell.interactionsOnsite, @SMatrix zeros(3, 3))
    push!(unitcell.interactionsField, @SVector zeros(3))
    return nothing
end

"""
    function addBasisSite!(unitcell::UnitCell{D}, position::NTuple{D,Float64}) where {D}
Convenience function that accepts position as an NTuple and converts it to an SVector internally.
"""
function addBasisSite!(unitcell::UnitCell{D}, position::NTuple{D,Float64}) where {D}
    addBasisSite!(unitcell, SVector(position))
end

"""
    function addBasisSite!(unitcell::UnitCell{D}, positions::Vector{SVector{D,Float64}}) where {D}
Convenience function that accepts multiple positions.
"""
function addBasisSite!(unitcell::UnitCell{D}, positions::Vector{SVector{D,Float64}}) where {D}
    for position in positions
        addBasisSite!(unitcell, position)
    end
    return nothing
end

"""
    function addBasisSite!(unitcell::UnitCell{D}, positions::Vector{NTuple{D,Float64}}) where {D}
Convenience function that accepts multiple positions and converts them into SVector internally.
"""
function addBasisSite!(unitcell::UnitCell{D}, positions::Vector{NTuple{D,Float64}}) where {D}
    for position in positions
        addBasisSite!(unitcell, SVector(position))
    end
    return nothing
end

"""
    function resetBasis!(unitcell::UnitCell{D}) where {D}
Resets the basis of a unitcell and all quantities that depend on the basis. 
This includes the interactions and interaction fields.
"""
function resetBasis!(unitcell::UnitCell{D}) where {D}
    #Also reset all quantities based on the basis
    unitcell.basis = Vector{SVector{D,Float64}}(undef, 0)
    unitcell.interactions = Vector{Tuple{Pair{Int,Int},NTuple{D,Int},SMatrix{3,3,Float64,9}}}(undef, 0)
    unitcell.interactionsOnsite = Vector{SMatrix{3,3,Float64,9}}(undef, 0)
    unitcell.interactionsField = Vector{SVector{3,Float64}}(undef, 0)
    return nothing
end

"""
--------------------------------------------------------------------------------
Extend Base
--------------------------------------------------------------------------------
"""

function Base.length(unitcell::UnitCell{D})::Int where {D}
    return length(unitcell.basis)
end

function Base.:show(io::IO, uc::UnitCell{D}) where {D}
    println(io, "$(dimension(uc))D Unitcell with $(length(uc.basis)) sites and $(length(uc.interactions)) interactions")
end

import Base: ==
function ==(uc1::UnitCell{D}, uc2::UnitCell{D}) where {D}
    return (
        uc1.primitive == uc2.primitive &&
        uc1.basis == uc2.basis &&
        uc1.interactions == uc2.interactions &&
        uc1.interactionsOnsite == uc2.interactionsOnsite &&
        uc1.interactionsField == uc2.interactionsField &&
        uc1.anisotropyFunction == uc2.anisotropyFunction &&
        uc1.anisotropyParameters == uc2.anisotropyParameters
    )
end