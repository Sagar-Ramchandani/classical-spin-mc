@testset "Monte Carlo Parameters" begin
    mcp = MonteCarloParameters(1.0, 10, 100)
    @test mcp.beta == 1.0
    @test mcp.thermalizationSweeps == 10
    @test mcp.measurementSweeps == 100

    """
    Test saving and loading Monte Carlo Parameters
    """

    h5open(tempname(), "w") do f
        @test writeMonteCarloParameters!(f, mcp) === nothing
        mcpLoaded = readMonteCarloParameters(f)
        @test typeof(mcp) == typeof(mcpLoaded)
        @test mcp == readMonteCarloParameters(f)
    end
end

@testset "Monte Carlo Statistics" begin
    mcs = MonteCarloStatistics()
    @test mcs.sweeps == 0

    """
    Test saving and loading Monte Carlo Statistics
    """

    h5open(tempname(), "w") do f
        @test writeMonteCarloStatistics!(f, mcs) === nothing
        @test mcs == readMonteCarloStatistics(f)
    end
end

@testset "Monte Carlo" begin
    mcp = MonteCarloParameters(1.0, 10, 100)
    mcp.beta == 1.0
    mcp.thermalizationSweeps == 10
    mcp.measurementSweeps == 100

    uc = UnitCell((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    addBasisSite!(uc, (0.5, 0.0, 0.0))
    addBasisSite!(uc, (0.0, 0.5, 0.0))
    addBasisSite!(uc, (0.0, 0.0, 0.5))
    addInteraction!(uc, 1 => 2, @SMatrix(ones(3, 3)), (0, 0, 0))
    addInteraction!(uc, 1, 3, @SMatrix(ones(3, 3)), (0, 0, 0))
    l = Lattice(uc, (3, 3, 3))

    mc = MonteCarlo(l, mcp)

    h5open(tempname(), "w") do f
        @test writeMonteCarlo!(f, mc) === nothing
        mcLoaded = readMonteCarlo(f)
        @test mcLoaded == mc
    end
end
