@testset "Lattice tests" begin

    #Generate a unitcell for testing purposes
    uc = UnitCell((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    addBasisSite!(uc, (0.5, 0.0, 0.0))
    addBasisSite!(uc, (0.0, 0.5, 0.0))
    addBasisSite!(uc, (0.0, 0.0, 0.5))
    addInteraction!(uc, 1 => 2, @SMatrix(ones(3, 3)), (0, 0, 0))
    addInteraction!(uc, 1, 3, @SMatrix(ones(3, 3)), (0, 0, 0))

    ll = Lattice(uc, (3, 3, 3))

    """
    Test saving and loading the lattice
    """

    h5open(tempname(), "w") do f
        @test writeLattice!(f, ll) === nothing
        @test ll == readLattice(f)
    end
end