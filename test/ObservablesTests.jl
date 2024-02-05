@testset "Observables" begin
    uc = UnitCell((1.0,))
    addBasisSite!(uc, (0.0,))
    l = Lattice(uc, (3,))
    obs = Observables(l)

    """
    Test saving and loading Observables 
    """

    h5open(tempname(), "w") do f
        @test writeObservables!(f, obs) === nothing
        obsLoaded = readObservables(f)
        @test typeof(obsLoaded) == typeof(obs)
        @test obs == obsLoaded
    end
end