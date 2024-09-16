Obsacc = 1e-14
include("UnitTests/UnitTests.jl")

using PMFRGSolve
using MPI
#@assert !isnothing(Base.get_extension(PMFRGSolve, :PMFRGSolveMPIExt)) "Perhaps you need `using MPI`?"
include("../ext/PMFRGSolveMPIExt/test/MPITest/mpi-tests.jl")

##
@testset verbose = true "PMFRGSolve tests" begin
    test_mpi()
    #testOneLoop(Obsacc)
    #testTwoLoop(Obsacc)
    #testParquet()
    #test_IO()
end