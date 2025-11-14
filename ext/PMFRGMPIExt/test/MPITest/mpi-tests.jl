using MPI
using Test

include("./test_partition.jl")
include("./test_best_partition_triangle.jl")
include("./test_best_partition_pairs.jl")

function test_mpi()
    @testset verbose = true "MPI tests" begin
        @testset verbose = true "Unit tests for MPI functionality" begin
            test_1D_partition()
            #test_best_partition_triangle()
            test_best_partition_pairs()
        end

        function run_mpi_script(script, n, testname)
            function print_header()
                linelength = 79
                println("="^linelength)
                title = " $testname "
                subtitle = " \"$(basename(script))\" "
                println("="^5 * title * "="^(linelength - 5 - length(title)))
                println("="^5 * subtitle * "="^(linelength - 5 - length(subtitle)))
                println("="^linelength)
            end

            print_header()
            @testset verbose = true "$testname" begin
                p = run(ignorestatus(`$(mpiexec())
                                      -n $n
                                      $(Base.julia_cmd())
                                      --project=$(Base.active_project())
                                      $script`))
                @test success(p)
            end
        end

        # DEBUG @testset verbose = true "MPI tests - external executables" begin
        # DEBUG     dir = dirname(@__FILE__)
        # DEBUG     run_mpi_script(
        # DEBUG         joinpath(dir, "..", "RegressionTests", "PMFRG.getXBubbleMPI.jl"),
        # DEBUG         2,
        # DEBUG         "Regression test - getXBubbleMPI",
        # DEBUG     )

        # DEBUG     run_mpi_script(
        # DEBUG         joinpath(dir, "test_chunk_communication.jl"),
        # DEBUG         4,
        # DEBUG         "Ibcast! communication example - test_chunk_communication.jl",
        # DEBUG     )

        # DEBUG     run_mpi_script(
        # DEBUG         joinpath(dir, "generate_data_example_mpi.jl"),
        # DEBUG         2,
        # DEBUG         "Generate Data Example",
        # DEBUG     )

        # DEBUG end
    end
end
