using StaticArrays
struct UnitCell{D}
    primitive::NTuple{D,SVector{D,Float64}}
    basis::Vector{SVector{D,Float64}}
    #interactions specified as (basis1=>basis2,offsetPrimitive,M)
    interactions::Vector{Tuple{Pair{Int,Int},NTuple{D,Int},SMatrix{3,3,Float64,9}}}
    interactionsOnsite::Vector{SMatrix{3,3,Float64,9}}
    interactionsField::Vector{SVector{3,Float64}}
end

#Constructor that accepts D number of SVectors as input to use as primitives. 
function UnitCell(primitive::Vararg{SVector{D,Float64},D}) where {D}
    return UnitCell{D}(primitive, Vector{SVector{D,Float64}}(undef, 0), Vector{Tuple{Pair{Int,Int},NTuple{D,Int},SMatrix{3,3,Float64,9}}}(undef, 0), Vector{SMatrix{3,3,Float64,9}}(undef, 0), Vector{SVector{3,Float64}}(undef, 0))
end

#Constructor with support for NTuples that get converted to SVectors.
function UnitCell(primitive::Vararg{NTuple{D,Float64},D}) where {D}
    return UnitCell(SVector.(primitive)...)
end

function Base.:show(io::IO, uc::UnitCell{D}) where {D}
    println(io, "$(2)D Unitcell with $(length(uc.basis)) sites and $(length(uc.interactions)) interactions")
end

"""
Adds an interaction between spin1 located at basis site `b1` of the given `unitcell` and spin2 at basis site `b2` in a unit cell that is offset by `offset` lattice vectors. 
The exchange energy is calculated as spin1'.M.spin2. 
"""
function addInteraction!(unitcell::UnitCell{D}, b::Pair{Int,Int}, M::SMatrix{3,3,Float64,9}, offset::NTuple{D,Int}=Tuple(zeros(Int, D))) where {D}
    if (b.first == b.second) && (offset == Tuple(zeros(Int, D)))
        @warn "Local Interaction detected, using setInteractionOnSite!() instead"
        setInteractionOnsite!(unitcell, b.first, M)
    end
    push!(unitcell.interactions, (b, offset, M))
    return nothing
end

function addInteraction!(unitcell::UnitCell{D}, b1::Int, b2::Int, M::SMatrix{3,3,Float64,9}, offset::NTuple{D,Int}=Tuple(zeros(Int, D))) where {D}
    addInteraction!(unitcell, b1 => b2, M, offset)
end

function addInteraction!(unitcell::UnitCell{D}, b1::Int, b2s::Vector{Int}, M::SMatrix{3,3,Float64,9}, offsets::Vector{NTuple{D,Int}}) where {D}
    for (b2, offset) in zip(b2s, offsets)
        addInteraction!(unitcell, b1 => b2, M, offset)
    end
    return nothing
end

"""
Sets the self-interaction for a spin located at basis site 'b' of the given 'unitcell'. 
The exchange energy is calculated as spin'.M.spin
"""
function setInteractionOnsite!(unitcell::UnitCell{D}, b::Int, M::SMatrix{3,3,Float64,9}) where {D}
    unitcell.interactionsOnsite[b] = M
    return nothing
end

"""
Sets the magnetic field for a spin at the basis site 'b'.
"""
function setField!(unitcell::UnitCell{D}, b::Int, B::SVector{3,Float64}) where {D}
    unitcell.interactionsField[b] = B
    return nothing
end

"""
Adds a new basis site to the given 'unitcell' at the given position with no self-interaction and no magnetic field.
"""
function addBasisSite!(unitcell::UnitCell{D}, position::SVector{D,Float64}) where {D}
    push!(unitcell.basis, position)
    push!(unitcell.interactionsOnsite, @SMatrix zeros(3, 3))
    push!(unitcell.interactionsField, @SVector zeros(3))
    return nothing
end

function addBasisSite!(unitcell::UnitCell{D}, position::NTuple{D,Float64}) where {D}
    addBasisSite!(unitcell, SVector(position))
end

function addBasisSites!(unitcell::UnitCell{D}, positions::Vector{SVector{D,Float64}}) where {D}
    for position in positions
        addBasisSite!(unitcell, position)
    end
    return nothing
end

function addBasisSites!(unitcell::UnitCell{D}, positions::Vector{NTuple{D,Float64}}) where {D}
    for position in positions
        addBasisSite!(unitcell, SVector(position))
    end
end

function resetBasis!(unitcell::UnitCell{D}) where {D}
    unitcell.basis = Vector{SVector{D,Float64}}(undef, 0)
    return nothing
end

function Base.length(unitcell::UnitCell{D})::Int where {D}
    return length(unitcell.basis)
end