# Examples

Here we showcase a few examples to demonstrate the 
use of this package.

## Defining a lattice spin model

As an example
```julia
using ClassicalSpinMC

#Create a unitcell from primitive lattice vectors a₁ and a₂.
a₁=(3/2, sqrt(3)/2)
a₂=(3/2, -sqrt(3)/2)
uc = UnitCell(a₁,a₂) 

#Add two basis sites to the unit cell at positions (0,0) and (1,0), respectively. 
addBasisSite!(uc, (0.0, 0.0))
addBasisSite!(uc, (1.0, 0.0)) 

#Add antiferromagnetic Heisenberg interactions between nearest neighbors.
M = [1.0 0.0 0.0;
    0.0 1.0 0.0;
    0.0 0.0 1.0] # Heisenberg interaction matrix

#Interaction of basis site 1 with site 2 in the same unit cell
addInteraction!(uc, 1=>2, M, (0, 0)) 

#Interaction of basis site 1 with site 2 in the unit cell shifted by (-1, 0) lattice vectors.
addInteraction!(uc, 1=>2, M, (-1, 0)) 

#Interaction of basis site b1 with site b2 in the unit cell shifted by (0, -1) lattice vectors.
addInteraction!(uc, 1=>2, M, (0, -1)) 

#Optionally apply a magnetic field B=(1,1,1) to basis site 1. 
#setField!(uc, 1, [1.0, 1.0, 1.0])

#Generate a lattice of 16*16 unit cells. 
L = (16, 16)
lattice = Lattice(uc, L)
```

