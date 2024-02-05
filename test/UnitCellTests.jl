@testset "UnitCell tests" begin

    """
    Checking creation of 1D, 2D and 3D UnitCell objects
    """

    @test typeof(UnitCell((1.0,))) == UnitCell{1}
    @test typeof(UnitCell((1.0, 0.0), (0.0, 1.0))) == UnitCell{2}
    @test typeof(UnitCell((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))) == UnitCell{3}

    """
    Further testing on 3D UnitCell object
    """

    uc = UnitCell((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    @test dimension(uc) == 3

    """
    Test adding a basis site
    """

    @test addBasisSite!(uc, (0.5, 0.0, 0.0)) === nothing
    @test addBasisSite!(uc, (0.0, 0.5, 0.0)) === nothing
    @test addBasisSite!(uc, (0.0, 0.0, 0.5)) === nothing

    """
    Test length of UnitCell fields
    """

    @test length(uc) == 3 #Tests for basis length
    @test length(uc.interactionsOnsite) == 3
    @test length(uc.interactionsField) == 3

    """
    Test for anisotropyFunction function
    """

    @test typeof(uc.anisotropyFunction) <: Function

    """
    Test for adding interactions
    """

    @test addInteraction!(uc, 1 => 2, @SMatrix(ones(3, 3)), (0, 0, 0)) === nothing
    @test addInteraction!(uc, 1, 3, @SMatrix(ones(3, 3)), (0, 0, 0)) === nothing
    @test length(uc.interactions) == 2

    """
    Test warning for onsite interactions added via addInteraction
    """

    @test_warn "Local Interaction detected, using setInteractionOnSite!() instead" addInteraction!(uc, 1, 1, @SMatrix(ones(3, 3)), (0, 0, 0))

    """
    Test setting onsite interaction for all sites
    """

    @test setInteractionOnsite!(uc, 1, 1.0 * @SMatrix(ones(3, 3))) === nothing
    @test setInteractionOnsite!(uc, 2, 2.0 * @SMatrix(ones(3, 3))) === nothing
    @test setInteractionOnsite!(uc, 3, 3.0 * @SMatrix(ones(3, 3))) === nothing

    @test uc.interactionsOnsite[1] == 1.0 * @SMatrix(ones(3, 3))
    @test uc.interactionsOnsite[2] == 2.0 * @SMatrix(ones(3, 3))
    @test uc.interactionsOnsite[3] == 3.0 * @SMatrix(ones(3, 3))

    """
    Test setting field interactions
    """

    @test setField!(uc, 1, 1.0 * @SVector(ones(3))) === nothing
    @test setField!(uc, 2, 2.0 * @SVector(ones(3))) === nothing
    @test setField!(uc, 3, 3.0 * @SVector(ones(3))) === nothing

    @test uc.interactionsField[1] == 1.0 * @SVector(ones(3))
    @test uc.interactionsField[2] == 2.0 * @SVector(ones(3))
    @test uc.interactionsField[3] == 3.0 * @SVector(ones(3))

    """
    Test saving and loading the unitcell
    """

    h5open(tempname(), "w") do f
        @test writeUnitcell!(f, uc) === nothing
        ucLoaded = readUnitcell(f)
        @test typeof(ucLoaded) == typeof(uc)
        @test uc == readUnitcell(f)
    end

    """
    Test resetting the basis
    """

    @test resetBasis!(uc) === nothing
    @test length(uc) == 0
    @test length(uc.interactions) == 0
    @test length(uc.interactionsOnsite) == 0
    @test length(uc.interactionsField) == 0
end