# This script can be run via
# ~/.julia/bin/mpiexecjl --project=./Playground
#                        -n 2
#                        julia <this file>
#

using Test
using HDF5
using MPI
using PMFRG
using SpinFRGLattices.SquareLattice

thisdir = dirname(@__FILE__)
non_mpi_regrtest_location = "../../../../test/RegressionTests"
include("../../../../test/RegressionTests/PMFRG.getXBubble.common.jl")


fname = joinpath(thisdir, non_mpi_regrtest_location, "PMFRG.getXBubble.data.h5")
h5file = h5open(fname, "r")

MPI.Init()

tests_ok = true
try
    @testset verbose = true "Tests for getXBubble!" begin
        ncases = read(h5file["Ncases"])
        @testset for i = 1:ncases
            (; X, State, Deriv, Lam) = h5deserialize(h5file, "arguments", i)
            X0 = X

            Par = generate_test_params()
            (; Buffs) = PMFRG.AllocateSetup(Par)
            Workspace = PMFRG.OneLoopWorkspace(State, Deriv, X0, Buffs, Par)

            PMFRG.getXBubble!(Workspace, Lam, UseMPI())

            (; X) = h5deserialize(h5file, "arguments_post", i)
            @test compare_arguments_post(X, X0)
        end
    end
catch e
    global tests_ok = false
    showerror(stdout, e)
end

MPI.Finalize()

exit((tests_ok ? 0 : 1))
