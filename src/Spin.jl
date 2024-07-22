"""
--------------------------------------------------------------------------------
Spin update functions 
--------------------------------------------------------------------------------
"""

"""
    function marsagliaSphereUpdate(rng=Random.GLOBAL_RNG)
Generates a random point on a unit sphere using the marsaglia method.
"""
function marsagliaSphereUpdate(rng = Random.GLOBAL_RNG)
    x1 = 0.0
    x2 = 0.0
    n = 0.0
    checkNormCondition = true

    while checkNormCondition
        x1 = 2 * rand(rng) - 1
        x2 = 2 * rand(rng) - 1
        n = x1^2 + x2^2
        if n < 1
            checkNormCondition = false
        end
    end

    return SVector(2 * x1 * sqrt(1 - n), 2 * x2 * sqrt(1 - n), 1 - 2 * n)
end

"""
    function sphericalUpdate(rng=Random.GLOBAL_RNG)
Generates a random point on a unit sphere using cylindrical co-ordinates.
"""
function sphericalUpdate(rng = Random.GLOBAL_RNG)
    phi = 2.0 * pi * rand(rng)
    z = 2.0 * rand(rng) - 1.0
    r = sqrt(1.0 - z * z)
    return SVector(r * cos(phi), r * sin(phi), z)
end

function conicalNorthPole(t::Float64, rng = Random.GLOBAL_RNG)
    phi = 2.0 * pi * rand(rng)
    zMinimal = cos(t)
    z = (1.0 - zMinimal) * rand(rng) + zMinimal
    return SVector(sqrt(1.0 - z^2) * cos(phi), sqrt(1.0 - z^2) * sin(phi), z)
end

"""
    function conicalUpdate(v::SVector{3,Float64}, t::Float64, rng=Random.GLOBAL_RNG)
Generates a random point on a unit sphere restricted in a cone defined as 
being centered around the vector `v` and an angle `t`.
For `t=pi` this is equivalent to any other sphere point picking method.
"""
function conicalUpdate(v::SVector{3, Float64}, t::Float64, rng = Random.GLOBAL_RNG)
    zAxis = SVector(0.0, 0.0, 1.0)
    c = dot(zAxis, v)
    s = sqrt(1 - c^2)
    rx, ry, rz = conicalNorthPole(t, rng)
    ux, uy, uz = normalize(cross(zAxis, v))
    rotatedVector = SVector(
        rx * (c + ux^2 * (1 - c)) + ry * (ux * uy * (1 - c) - (uz * s)) +
        rz * (ux * uz * (1 - c) + (uy * s)),
        rx * (ux * uy * (1 - c) + (uz * s)) + ry * (c + (uy^2) * (1 - c)) +
        rz * (uy * uz * (1 - c) - (ux * s)),
        rx * (ux * uz * (1 - c) - (uy * s)) + ry * (uy * uz * (1 - c) + (ux * s)) +
        rz * (c + (uz^2) * (1 - c))
    )
    return rotatedVector
end

"""
--------------------------------------------------------------------------------
Energy functions 
--------------------------------------------------------------------------------
"""

"""
    function exchangeEnergy(s1::SVector{3,Float64}, M::SMatrix{3,3,Float64,9}, s2::SVector{3,Float64})
Calculates the energy contribution due to exchange terms.
"""
function exchangeEnergy(
        s1::SVector{3, Float64}, M::SMatrix{3, 3, Float64, 9}, s2::SVector{3, Float64})
    return (s1' * M) * s2
end

"""
    function getAnisotropy(::typeof(identity), lattice::Lattice{D,N})::Float64 where {D,N}
Default method for returning zero in case of no-anisotropy for the full lattice.
"""
function getAnisotropy(::typeof(identity), lattice::Lattice{D, N})::Float64 where {D, N}
    return zero(Float64)
end

"""
    function getAnisotropy(::typeof(identity), spin::SVector{3,Float64}, parameters::NTuple{P,Float64})::Float64 where {P}
Default method for returning zero in case of no-anisotropy for a particular spin.
"""
function getAnisotropy(::typeof(identity), spin::SVector{3, Float64},
        parameters::NTuple{P, Float64})::Float64 where {P}
    return zero(Float64)
end

"""
    function getAnisotropy(func::F, spin::SVector{3,Float64}, parameters::NTuple{P,Float64})::Float64 where {F<:Function,P}
Calculates the anisotropy energy term for a particular spin.
"""
function getAnisotropy(func::F, spin::SVector{3, Float64},
        parameters::NTuple{P, Float64})::Float64 where {F <: Function, P}
    return func(spin, parameters)
end

"""
    function getAnisotropy(func::F, lattice::Lattice{D,N})::Float64 where {D,N,F<:Function}
Calculates the anisotropy energy term for the full lattice.
"""
function getAnisotropy(func::F, lattice::Lattice{D, N})::Float64 where {D, N, F <: Function}
    energy = zero(Float64)
    for spin in lattice.spins
        energy += getAnisotropy(func, spin, lattice.anisotropyParameteres)
    end
    return energy
end

"""
    function getEnergy(lattice::Lattice{D,N})::Float64 where {D,N}
Calculates the energy of the current spin configuration of the lattice 
based on all terms defined such as exchange, onsite, field and ansiotropy terms.
"""
function getEnergy(lattice::Lattice{D, N})::Float64 where {D, N}
    energy = zero(Float64)
    for site in 1:length(lattice)
        s0 = getSpin(lattice, site)

        #two-spin interactions
        interactionSites = getInteractionSites(lattice, site)
        interactionMatrices = getInteractionMatrices(lattice, site)
        for (i, j) in zip(interactionMatrices, interactionSites)
            if site > j
                energy += exchangeEnergy(s0, i, getSpin(lattice, j))
            end
        end

        #onsite interaction
        energy += exchangeEnergy(s0, getInteractionOnsite(lattice, site), s0)

        #field interaction
        energy += dot(s0, getInteractionField(lattice, site))

        #ansiotropy
        energy += getAnisotropy(lattice.anisotropyFunction, lattice)
    end

    return energy
end

"""
    function getEnergyDifference(lattice::Lattice{D,N}, site::Int, newState::SVector{3,Float64})::Float64 where {D,N}
Performant method for calculating the energy difference when updating a `site` into a `newState`
due to all terms such as exchange, onsite, field and ansiotropy terms.
"""
function getEnergyDifference(lattice::Lattice{D, N}, site::Int,
        newState::SVector{3, Float64})::Float64 where {D, N}
    ΔE = 0.0
    oldState = getSpin(lattice, site)
    ΔS = newState - oldState

    #two-spin interactions
    interactionSites = getInteractionSites(lattice, site)
    interactionMatrices = getInteractionMatrices(lattice, site)
    for (i, j) in zip(interactionMatrices, interactionSites)
        ΔE += exchangeEnergy(ΔS, i, getSpin(lattice, j))
    end

    #onsite interaction
    interactionOnsite = getInteractionOnsite(lattice, site)
    ΔE += exchangeEnergy(newState, interactionOnsite, newState) -
          exchangeEnergy(oldState, interactionOnsite, oldState)

    #field interaction
    ΔE += dot(ΔS, getInteractionField(lattice, site))

    #ansiotropy
    ΔE += getAnisotropy(
        lattice.anisotropyFunction, newState, lattice.anisotropyParameteres) -
          getAnisotropy(lattice.anisotropyFunction, oldState, lattice.anisotropyParameteres)

    return ΔE
end

"""
--------------------------------------------------------------------------------
Observables functions 
--------------------------------------------------------------------------------
"""

"""
    function getMagnetization(lattice::Lattice{D,N}) where {D,N}
Measures the Magnetization of the lattice spin configuration.
"""
function getMagnetization(lattice::Lattice{D, N}) where {D, N}
    return mean(lattice.spins)
end

"""
    function getMagnetizationPerSite(lattice::Lattice{D,N}) where {D,N}
Measures the Magnetization for each type of basis site in the unitcell seperately.
"""
function getMagnetizationPerSite(lattice::Lattice{D, N}) where {D, N}
    nTypes = length(lattice.unitcell)
    mSites = zeros(3, nTypes)

    for ii in 1:nTypes
        total = (0.0, 0.0, 0.0)
        for jj in range(ii, lattice.length, step = nTypes)
            total = total .+ getSpin(lattice, jj)
        end
        mSites[:, ii] .= total
    end
    mSites ./= lattice.length / nTypes
    return mSites
end

"""
    function calcTriangles(lattice::Lattice{D,N}) where {D,N}
Returns a list of all unique triangles found in the lattice
consisting of the sites comprising the triangle in ascending order.
"""
function calcTriangles(lattice::Lattice{D, N}) where {D, N}
    siteConnections = lattice.interactionSites
    nTypes = length(lattice.unitcell.basis)
    triangles = Vector{Vector{Int64}}()
    for ii in 1:(lattice.length)
        allConnections = siteConnections[ii]
        secondLevelConnections = [siteConnections[jj] for jj in allConnections]
        for jj in allConnections
            currentIndex = jj
            for kk in eachindex(secondLevelConnections)
                if currentIndex in secondLevelConnections[kk]
                    tri = [ii, currentIndex, allConnections[kk]]
                    triType = [ii % nTypes for ii in tri]
                    perm = sortperm(triType)
                    push!(triangles, tri[perm])
                end
            end
        end
    end
    result = (unique(triangles))
    return result
end

"""
    function getChirality(lattice::Lattice{D,N}, siteList::Vector{Vector{Int64}}) where {D,N}
Given a siteList of triangles in the lattice, computes the chirality of the lattice,
where the chirality is defined as the vector triple product of spins on a triangle.
"""
function getChirality(lattice::Lattice{D, N}, siteList::Vector{Vector{Int64}}) where {D, N}
    chi = 0
    for ii in siteList
        chi += dot(getSpin(lattice, ii[1]),
            cross(getSpin(lattice, ii[2]), getSpin(lattice, ii[3])))
    end
    chi /= length(siteList)
    return chi
end

"""
    function getCorrelation(lattice::Lattice{D,N}) where {D,N}
Measures the two-point correlation for each pair of sites in the lattice.
"""
function getCorrelation(lattice::Lattice{D, N}) where {D, N}
    corr = zeros(length(lattice), length(lattice.unitcell.basis))
    for i in 1:length(lattice.unitcell.basis)
        s0 = getSpin(lattice, i)
        for j in 1:length(lattice)
            corr[j, i] = dot(s0, getSpin(lattice, j))
        end
    end
    return corr
end

"""
--------------------------------------------------------------------------------
Overrelaxation functions 
--------------------------------------------------------------------------------
"""

"""
    function microcanonicalRotation(lattice::Lattice{D,N}, site::Int) where {D,N}
Performs an Overrelaxation step on the lattice in a deterministic way by 
using a π rotation.

!!! note "Assumptions on interactions and ansiotropy"
    This assumes that there are no anisotropy contributions 
    and no self interactions for the lattice object. 
    For performance reasons, these checks are performed seperately.
"""
function microcanonicalRotation(lattice::Lattice{D, N}, site::Int) where {D, N}
    #compute rotation axis
    axis = SVector(0.0, 0.0, 0.0)
    for (interactionSite, interactionMatrix) in zip(
        lattice.interactionSites[site], lattice.interactionMatrices[site])
        axis = axis .+ interactionMatrix * getSpin(lattice, interactionSite)
    end
    axis = axis .+ getInteractionField(lattice, site)
    axis = normalize(axis)

    #pi-rotation of current spin around axis
    return piRotation(getSpin(lattice, site), axis)
end

"""
    function microcanonicalRotationRandom(lattice::Lattice{D,N}, site::Int, rng=Random.GLOBAL_RNG) where {D,N}
Performs an Overrelaxation step on the lattice in a non-deterministic way by 
using a random rotation.

!!! note "Assumptions on interactions and anisotropy"
    This assumes that there are no anisotropy contributions 
    and no self interactions for the lattice object. 
    For performance reasons, these checks are performed seperately.
"""
function microcanonicalRotationRandom(
        lattice::Lattice{D, N}, site::Int, rng = Random.GLOBAL_RNG) where {D, N}
    #compute rotation axis
    axis = zero(eltype(lattice.spins))
    for (interactionSite, interactionMatrix) in zip(
        lattice.interactionSites[site], lattice.interactionMatrices[site])
        axis = axis .+ interactionMatrix * getSpin(lattice, interactionSite)
    end
    axis = axis .+ getInteractionField(lattice, site)
    axis = normalize(axis)

    #random-rotation of current spin around axis
    return rotateAboutAxis(getSpin(lattice, site), axis, 2π * rand(rng))
end

"""
    function rotateAboutAxis(spin::SVector{3,Float64}, axis::SVector{3,Float64}, angle::Float64)
Rotates a vector about an axis by a particular angle. 

!!! note "Assumptions on normalization"
    This assumes that the vector and axis are normalized 
    for performance reasons.
"""
function rotateAboutAxis(
        spin::SVector{3, Float64}, axis::SVector{3, Float64}, angle::Float64)
    spinParallel = dot(spin, axis) .* axis
    spinPerp = spin .- spinParallel
    w = cross(axis, spinPerp)
    x1 = cos(angle) ./ norm(spinPerp)
    x2 = sin(angle) ./ norm(w)
    spinPerpRotated = norm(spinPerp) .* ((x1 .* spinPerp) .+ (x2 .* w))
    return (spinPerpRotated .+ spinParallel)
end

"""
    function piRotation(spin::SVector{3,Float64}, axis::SVector{3,Float64})
Special case of `rotateAboutAxis` for a rotation by an angle π.

!!! note "Assumptions on normalization"
    This assumes that the vector and axis are normalized 
    for performance reasons.
"""
function piRotation(spin::SVector{3, Float64}, axis::SVector{3, Float64})
    return SVector(
        2 * (axis[1] * axis[2] * spin[2] + axis[1] * axis[3] * spin[3] +
         axis[1] * axis[1] * spin[1]) - spin[1],
        2 * (axis[1] * axis[2] * spin[1] + axis[2] * axis[3] * spin[3] +
         axis[2] * axis[2] * spin[2]) - spin[2],
        2 * (axis[1] * axis[3] * spin[1] + axis[2] * axis[3] * spin[2] +
         axis[3] * axis[3] * spin[3]) - spin[3]
    )
end
