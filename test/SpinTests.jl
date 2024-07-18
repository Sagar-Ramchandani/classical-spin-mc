@testset "Spin tests" begin
    @testset "Spherical Distribution of Spins:sphericalUpdate" begin
        N = 10_000_000
        total = @SVector(zeros(3))
        for _ in 1:N
            total += sphericalUpdate()
        end
        @test all(isapprox.(total ./ N, zeros(3), atol = 1.0e-3))
    end
    @testset "Spherical Distribution of Spins:marsagliaSphereUpdate" begin
        N = 10_000_000
        total = @SVector(zeros(3))
        for _ in 1:N
            total += marsagliaSphereUpdate()
        end
        @test all(isapprox.(total ./ N, zeros(3), atol = 1.0e-3))
    end
    @testset "Spherical Distribution of Spins:conicalUpdate" begin
        N = 10_000_000
        total = @SVector(zeros(3))
        for _ in 1:N
            total += conicalUpdate(@SVector(ones(3)), 1.0π)
        end
        @test all(isapprox.(total ./ N, zeros(3), atol = 1.0e-3))
    end
    @testset "Magnitude of Spins:sphericalUpdate" begin
        N = 10_000_000
        testPass = true
        for _ in 1:N
            if !(isapprox(norm(sphericalUpdate()), 1.0))
                testPass = false
                break
            end
        end
        @test testPass
    end
    @testset "Magnitude of Spins:marsagliaSphereUpdate" begin
        N = 10_000_000
        testPass = true
        for _ in 1:N
            if !(isapprox(norm(marsagliaSphereUpdate()), 1.0))
                testPass = false
                break
            end
        end
        @test testPass
    end
    @testset "Magnitude of Spins:conicalUpdate" begin
        N = 10_000_000
        testPass = true
        for _ in 1:N
            if !(isapprox(norm(conicalUpdate(@SVector(ones(3)), rand())), 1.0))
                testPass = false
            end
        end
        @test testPass
    end

    """
    Testing exchange energy
    """

    M = SMatrix{3, 3, Float64, 9}([1.1 2.2 3.3; 4.4 5.5 6.6; 7.7 8.8 9.9])
    s1 = SVector(1.0, 2.0, 3.0)
    s2 = SVector(4.0, 5.0, 6.0)
    @test ClassicalSpinMC.exchangeEnergy(s1, M, s2) ≈ 607.2

    """
    Create a lattice for further testing
    """

    a1 = (3 / 2, sqrt(3) / 2)
    a2 = (3 / 2, -sqrt(3) / 2)
    basis = [(0.0, 0.0), (1.0, 0.0)]
    uc = UnitCell(a1, a2)
    b1 = addBasisSite!(uc, basis[1])
    b2 = addBasisSite!(uc, basis[2])
    addInteraction!(uc, 1, 2,
        SMatrix{3, 3, Float64, 9}([-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]), (0, 0))
    addInteraction!(uc, 1, 2,
        SMatrix{3, 3, Float64, 9}([-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]), (0, -1))
    addInteraction!(uc, 1, 2,
        SMatrix{3, 3, Float64, 9}([-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]), (-1, 0))
    lattice = Lattice(uc, (4, 4))

    """
    Setting spins to predetermined configurations
    """

    for i in 1:length(lattice)
        setSpin!(lattice, i, (0.0, 0.0, 1.0))
    end

    """
    Test lattice parameters
    """

    @test getEnergy(lattice) ≈ -48.0
    @test ClassicalSpinMC.getEnergyDifference(lattice, 1, SVector(0.0, 0.0, -1.0)) ≈ 6.0
    setSpin!(lattice, 1, (0.0, 0.0, -1.0))
    @test getEnergy(lattice) ≈ -42.0

    @test getMagnetization(lattice) ≈ SVector(0.0, 0.0, 30.0 / 32.0)

    @test all(getCorrelation(lattice)[:, 1] .≈ [(i == 1 ? 1.0 : -1.0) for i in 1:32])
    @test all(getCorrelation(lattice)[:, 2] .≈ [(i == 1 ? -1.0 : 1.0) for i in 1:32])
end
