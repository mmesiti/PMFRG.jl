using MPI
using Test

include("./test_partition.jl")
include("./test_best_partition_pairs.jl")
function test_mpi_core()
    @testset verbose = true "MPI tests" begin
        @testset verbose = true "Unit tests for MPI functionality" begin
            test_1D_partition()
            test_best_partition_pairs()
        end


        function run_mpi_script(script, n, testname)
            linelength = 79
            function print_header()
                println("="^linelength)
                title = " $testname "
                subtitle = " \"$(basename(script))\" "
                function print_line(string)
                    println("="^5 * string * "="^(linelength - 5 - length(string)))
                end

                print_line(title)
                print_line(subtitle)
                print_line(" Output of first rank only ")
                println("="^linelength)
            end
            print_footer() = println("="^linelength)

            function create_mpi_shell_wrapper()
                # From: https://www.open-mpi.org/community/lists/users/2012/02/18362.php
                mpi_shell_wrapper_text = """
#!/bin/bash
ARGS=\$@
if [[ -v OMPI_COMM_WORLD_RANK ]] # OpenMPI Case
then
   MY_MPI_RANK=\$OMPI_COMM_WORLD_RANK
elif [[ -v PMI_RANK ]]           # MPICH Case, probably
then
   MY_MPI_RANK=\$PMI_RANK
fi
if [[ \$MY_MPI_RANK == 0 ]]
then
  \$ARGS
else
  \$ARGS 1>/dev/null 2>/dev/null
fi
"""
                path = "./shell-wrapper-text.sh"
                open(path, "w") do io
                    write(io, mpi_shell_wrapper_text)
                end
                run(`chmod +x $path`)
                path
            end


            print_header()
            wrapper_path = create_mpi_shell_wrapper()
            @testset verbose = true "$testname" begin
                p = run(ignorestatus(`$(mpiexec())
                                      -n $n
                                      $wrapper_path
                                      $(Base.julia_cmd())
                                      --project=$(Base.active_project())
                                      $script`))
                @test success(p)
            end
            run(`rm $wrapper_path`)
            print_footer()
        end

        << << << < HEAD:PMFRGCore.jl/ext/PMFRGCoreMPIExt/test/MPITest/mpi-tests.jl
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
        # DEBUG        joinpath(dir, "test_setup_functions.jl"),
        # DEBUG        2,
        # DEBUG        "Setup and compatibility tests - MPI version",
        # DEBUG     )
        # DEBUG     run_mpi_script(
        # DEBUG         joinpath(dir, "generate_data_example_mpi.jl"),
        # DEBUG         2,
        # DEBUG         "Generate Data Example",
        # DEBUG     )
        # DEBUG end
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
