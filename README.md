# ClassicalSpinMC.jl

*Markov chain Monte Carlo for 
classical lattice spin models*

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://sagar-ramchandani.github.io/classical-spin-mc/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sagar-ramchandani.github.io/classical-spin-mc/dev/)

This package aims to provide a flexible and performant implementation of the MCMC algorithm for classical lattice spin models.
This project is a fork of [SpinMC.jl](https://github.com/fbuessen/SpinMC.jl)

## Package features at a glance
- Simulated annealing.
- Replica exchange / Parallel tempering.
- Arbitrary D-dimensional lattices.
- Microcanonical / Over relaxation sweeps.
- Conical spin-flip updates.
- Supports user-defined update methods and observables
- Data storage in HDF5 format.

## Brief overview

The package can simulate spin models described by the Hamiltonian,\
$H=\sum_{i,j} S_i \cdot M_{ij} \cdot S_j + \sum_{i} S_i \cdot B_j +
A(S_i)$\
Here $S_i$ are classical $O(3)$ spins, with interactions $M$, 
magnetic field $B$ and an arbitrary anisotropy function $A$.

By default, the package has the following methods for updating spins.
- Marsaglia method.
- Cylindrical co-ordinates.
- Conical updates.

Conical updates sample around the current position of the spin on the unit sphere to maintain an acceptance rate of 50% thereby allowing for significantly reduced auto-correlation time for low temperature measurements. 

The statistics for the measurements are provided using the
[BinningAnalysis.jl](https://github.com/carstenbauer/BinningAnalysis.jl) package. By default, the package measures energy, specific heat, 
magnetization and spin-spin correlation. This behavior can be extended, please see ... for this.

For better tempering, the package implements simulated annealing 
and parallel tempering. The latter of which uses [Distributed.jl](https://github.com/JuliaLang/Distributed.jl) allowing the sampling of large number of temperature points across multiple machines such 
as on a cluster.
