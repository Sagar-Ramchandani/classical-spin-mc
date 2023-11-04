using StaticArrays
mutable struct Lattice{D,N}
    size::NTuple{D,Int} #linear dimension (D) of the lattice in number of unit cells
    length::Int #Number of sites in the lattice N_sites
    unitcell::UnitCell{D}
    sitePositions::Vector{SVector{D,Float64}}

    spins::Vector{SVector{3,Float64}} #Vector of spins

    interactionSites::Vector{NTuple{N,Int}} #list of length N_sites, for every site contains all interacting sites
    interactionMatrices::Vector{NTuple{N,SMatrix{3,3,Float64,9}}} #list of length N_sites, for every site contains all interaction matrices
    interactionOnsite::Vector{SMatrix{3,3,Float64,9}} #list of length N_sites, for every site contains the local onsite interaction matrix
    interactionField::Vector{SVector{3,Float64}} #list of length N_sites, for every site contains the local field

    #Generic fallback method, possibly remove
    Lattice(D, N) = new{D,N}()
end

#Add error check for a unitcell with no sites.
function Lattice(uc::UnitCell{D}, L::NTuple{D,Int}) where {D}
    #parse interactions
    ##For every basis site b, generate list of sites which b interacts with and store the corresponding interaction sites and matrices. 
    ##Interaction sites are specified by the target site's basis id, b_target, and the offset in units of primitive lattice vectors. 
    ##If b has multiple interactions defined with the same target site, eliminate those duplicates by summing up the interaction matrices. 
    interactionTargetSites = [Vector{Tuple{Int,NTuple{D,Int},SMatrix{3,3,Float64,9}}}(undef, 0) for i in 1:length(uc.basis)] #tuples of (b_target, offset, M)
    for x in uc.interactions
        b, offset, M = x
        b1, b2 = b
        if b1 == b2 && offset .% L == Tuple(zeros(D))
            #Possibly automatically add using setInteractionOnSite with a warning?
            error("Interaction cannot be local. Use setInteractionOnsite!() instead.")
        end

        #locate existing coupling to target site and add interaction matrix
        for i in 1:length(interactionTargetSites[b1])
            if interactionTargetSites[b1][i][1] == b2 && interactionTargetSites[b1][i][2] == offset
                interactionTargetSites[b1][i] = (interactionTargetSites[b1][i][1], interactionTargetSites[b1][i][2], interactionTargetSites[b1][i][3] + M)
                @goto endb1
            end
        end
        #if coupling does not exist yet, push new entry
        push!(interactionTargetSites[b1], (b2, offset, M))
        @label endb1

        #locate existing coupling from target site and add interaction matrix
        for i in 1:length(interactionTargetSites[b2])
            if interactionTargetSites[b2][i][1] == b1 && interactionTargetSites[b2][i][2] == (x -> -x).(offset)
                interactionTargetSites[b2][i] = (interactionTargetSites[b2][i][1], interactionTargetSites[b2][i][2], interactionTargetSites[b2][i][3] + transpose(M))
                @goto endb2
            end
        end
        #if coupling does not exist yet, push new entry
        push!(interactionTargetSites[b2], (b1, (x -> -x).(offset), transpose(M)))
        @label endb2
    end
    Ninteractions = findmax([length(interactionTargetSites[i]) for i in 1:length(uc.basis)])[1]

    #create lattice struct
    lattice = Lattice(D, Ninteractions)
    lattice.size = L
    lattice.length = prod(L) * length(uc.basis)
    lattice.unitcell = uc

    #generate linear representation of lattice sites to assign integer site IDs
    ##Enumeration sequence is (a1, a2, ..., b) in row-major fashion
    sites = Vector{NTuple{D + 1,Int}}(undef, lattice.length)
    function nextSite(site)
        next = collect(site)
        next[D+1] += 1
        if next[D+1] > length(uc.basis)
            next[D+1] = 1
            next[D] += 1
        end
        for d in reverse(1:D)
            if next[d] >= L[d]
                next[d] = 0
                d - 1 > 0 && (next[d-1] += 1)
            end
        end
        return Tuple(next)
    end
    sites[1] = tuple(zeros(Int, D)..., 1)
    for i in 2:length(sites)
        sites[i] = nextSite(sites[i-1])
    end

    #init site positions
    lattice.sitePositions = Vector{NTuple{D,Float64}}(undef, length(sites))
    for i in 1:length(sites)
        site = sites[i]
        lattice.sitePositions[i] = .+([uc.primitive[j] .* site[j] for j in 1:D]...) .+ uc.basis[site[end]]
    end

    #init spins 
    lattice.spins = Vector{SVector{3,Float64}}(undef, lattice.length)

    #write interactions to lattice
    lattice.interactionSites = repeat([NTuple{Ninteractions,Int}(ones(Int, Ninteractions))], lattice.length)
    lattice.interactionMatrices = repeat([NTuple{Ninteractions,SMatrix{3,3,Float64,9}}(repeat([@SMatrix zeros(3, 3)], Ninteractions))], lattice.length)
    lattice.interactionOnsite = repeat([@SMatrix zeros(3, 3)], lattice.length)
    lattice.interactionField = repeat([@SVector zeros(3)], lattice.length)

    function applyPBC(n, L)
        while n < 0
            n += L
        end
        while n >= L
            n -= L
        end
        return n
    end
    function siteIndexFromParametrization(site)
        return findfirst(isequal(site), sites)
    end

    for i in 1:length(sites)
        site = sites[i]
        b = site[end]

        #onsite interaction
        lattice.interactionOnsite[i] = uc.interactionsOnsite[b]

        #field interaction
        lattice.interactionField[i] = uc.interactionsField[b]

        #two-spin interactions
        interactionSites = repeat([i], Ninteractions)
        interactionMatrices = repeat([@SMatrix zeros(3, 3)], Ninteractions)
        for j in 1:Ninteractions
            if j <= length(interactionTargetSites[b])
                b2, offset, M = interactionTargetSites[b][j]

                primitiveTarget = [applyPBC(site[k] + offset[k], L[k]) for k in 1:D]
                targetSite = tuple(primitiveTarget..., b2)

                interactionSites[j] = siteIndexFromParametrization(targetSite)
                interactionMatrices[j] = InteractionMatrix(M)
            end
        end
        lattice.interactionSites[i] = NTuple{Ninteractions,Int}(interactionSites)
        lattice.interactionMatrices[i] = NTuple{Ninteractions,SMatrix{3,3,Float64,9}}(interactionMatrices)
    end
    return lattice
end

function Base.size(lattice::Lattice{D,N}) where {D,N}
    return lattice.size
end

function Base.length(lattice::Lattice{D,N}) where {D,N}
    return lattice.length
end

function Base.:show(io::IO, lattice::Lattice{D,N}) where {D,N}
    println(io, "$(D)D Lattice with $(size(lattice)) unitcells and $(N) interactions per site")
end

function getSpin(lattice::Lattice{D,N}, site::Int) where {D,N}
    return lattice.spins[site]
end

function setSpin!(lattice::Lattice{D,N}, site::Int, newState::SVector{3,Float64}) where {D,N}
    lattice.spins[site] = newState
end

function getSitePosition(lattice::Lattice{D,N}, site::Int)::SVector{D,Float64} where {D,N}
    return lattice.sitePositions[site]
end

function getInteractionSites(lattice::Lattice{D,N}, site::Int)::NTuple{N,Int} where {D,N}
    return lattice.interactionSites[site]
end

function getInteractionMatrices(lattice::Lattice{D,N}, site::Int)::NTuple{N,SMatrix{3,3,Float64,9}} where {D,N}
    return lattice.interactionMatrices[site]
end

function getInteractionOnsite(lattice::Lattice{D,N}, site::Int)::SMatrix{3,3,Float64,9} where {D,N}
    return lattice.interactionOnsite[site]
end

function getInteractionField(lattice::Lattice{D,N}, site::Int)::SVector{3,Float64} where {D,N}
    return lattice.interactionField[site]
end