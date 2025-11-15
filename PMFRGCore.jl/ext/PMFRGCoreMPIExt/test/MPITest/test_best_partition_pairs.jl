using Test
include("../../mpi/best_partition_pairs.jl")

using .BestPartitionPairs:
    _get_ranges_t, _all_ns_x_nt_factorizations, get_all_ranges_st, get_imbalance_from_ranges


function test_best_partition_pairs()
    @testset verbose = true "Tests for best it pairs partition" begin
        @testset verbose = true "Simple manual cases" begin
            #                    N,nranks,rank
            @test _get_ranges_t(10, 2, 0) in [(1:3, 9:10), (1:2, 8:10)]

            @test _get_ranges_t(10, 2, 1) in [(4:8,), (3:7,)]

            @test _get_ranges_t(11, 2, 0) in [(1:2, 9:11), (1:3, 10:11), (1:3, 9:11)]

            @test _get_ranges_t(11, 2, 1) in [(3:8,), (4:9,), (4:8,)]

            @test _get_ranges_t(12, 3, 0) == (1:2, 11:12)
            @test _get_ranges_t(12, 3, 1) == (3:4, 9:10)
            @test _get_ranges_t(12, 3, 2) == (5:8,)

        end
        @testset "ranges cover all sites" begin
            for N = 2:30, nranks = 1:min(5, div(N, 2))
                coverage = zeros(Int, N)
                all_ranges_for_all_ranks =
                    [_get_ranges_t(N, nranks, irank) for irank = 0:(nranks-1)]
                for rank = 1:nranks
                    for itrange in all_ranges_for_all_ranks[rank]
                        for it in itrange
                            coverage[it] += 1
                        end
                    end
                end
                for it = 1:N
                    if coverage[it] != 1
                        println("$nranks ransk, $N freqs")
                        println(all_ranges_for_all_ranks)

                    end
                    @test coverage[it] == 1
                end
            end
        end

        @testset "ranges do not go out of 1:N range" begin
            for N = 2:30, nranks = 1:min(5, div(N, 2))
                all_ranges_for_all_ranks =
                    [_get_ranges_t(N, nranks, irank) for irank = 0:(nranks-1)]
                for rank = 1:nranks
                    for itrange in all_ranges_for_all_ranks[rank]
                        @test 1 <= itrange.start <= N
                        @test 1 <= itrange.stop <= N
                    end
                end
            end
        end

        @testset "_get_ranges_t requires N >= 2*nranks_t" begin
            # Test that _get_ranges_t throws TooManyRanksError when N < 2*nranks_t
            @test_throws BestPartitionPairs.TooManyRanksError _get_ranges_t(3, 2, 0)  # N < 2*nranks_t
            @test _get_ranges_t(4, 2, 0) isa Tuple  # N = 2*nranks_t (boundary case)
            @test _get_ranges_t(10, 2, 0) isa Tuple # N > 2*nranks_t
        end

        @testset "get_all_ranges_st covers all (is, it) pairs exactly once" begin
            for N = 10:13, nranks = 2:min(5, N)
                all_ranges = get_all_ranges_st(N, nranks)

                # Check correct number of ranks
                @test length(all_ranges) == nranks

                # Check all (is, it) pairs covered exactly once
                coverage = zeros(Int, N, N)
                for (isrange, itranges) in all_ranges
                    for is in isrange
                        for itrange in itranges
                            for it in itrange
                                coverage[is, it] += 1
                            end
                        end
                    end
                end
                @test all(coverage .== 1)
            end
        end

        @testset "N=50 with nranks=2,4,8 - coverage and load balancing" begin
            N = 50
            for nranks in [2, 4, 8]
                @testset "N=$N, nranks=$nranks" begin
                    all_ranges = get_all_ranges_st(N, nranks)

                    # Check correct number of ranks
                    @test length(all_ranges) == nranks

                    # Check all (is, it) pairs covered exactly once
                    coverage = zeros(Int, N, N)
                    for (isrange, itranges) in all_ranges
                        for is in isrange
                            for itrange in itranges
                                for it in itrange
                                    coverage[is, it] += 1
                                end
                            end
                        end
                    end
                    @test all(coverage .== 1)

                    # Check load balancing metrics
                    imbalance_X, imbalance_XTilde = get_imbalance_from_ranges(all_ranges, N)

                    # Verify reasonable load balancing (imbalance should be < 50%)
                    @test imbalance_X < 0.07
                    @test imbalance_XTilde < 0.07
                end
            end
        end
    end
end
