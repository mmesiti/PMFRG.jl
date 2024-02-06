Obsacc = 1e-14
include("UnitTests/UnitTests.jl")

using MPI, PencilArrays
@assert !isnothing(Base.get_extension(PMFRG, :PMFRGMPIExt)) "Perhaps you need `using MPI, PencilArrays`?"
include("../ext/PMFRGMPIExt/test/MPITest/mpi-tests.jl")
include("RegressionTests/PMFRG.getXBubble.jl")

@testset verbose = true "PMFRG tests" begin
    testStateUnpacking()
    testOneLoop(Obsacc)
    testTwoLoop(Obsacc)
    testParquet()
    test_IO()
    test_mpi()
    test_getXBubble()
end
