using Test
include("../../mpi/best_partition_pairs.jl")

using .BestPartitionPairs:
    _get_ranges_t,
    _get_number_of_sites_eo,
    _all_ns_x_ntu_factorizations,
    get_all_ranges_stu,
    get_imbalance


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
            for N = 1:30, nranks = 1:min(5, N)
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

        #@testset "ranges do not go out of 1:N range" begin
        #    for N = 1:10, nranks = 1:min(5, N)
        #        all_ranges = [_get_ranges_tu(N, nranks, irank) for irank = 0:(nranks-1)]
        #        for (itrange, iurange) in all_ranges
        #            @test 1 <= itrange.start <= N
        #            @test 1 <= itrange.stop <= N
        #            @test 1 <= iurange.start <= N
        #            @test 1 <= iurange.stop <= N
        #        end
        #    end
        #end




        #@testset "get all nranks_s x nranks_tu factorizations" begin
        #    all_pairs = _all_ns_x_ntu_factorizations(50)
        #    @test (1, 50) in all_pairs
        #    @test (2, 25) in all_pairs
        #    for pair in all_pairs
        #        @test typeof(pair) <: Tuple{Int,Int}
        #    end
        #    for i = 1:50
        #        if 50 % i != 0
        #            for pair in all_pairs
        #                @test pair[1] != i
        #            end
        #        end
        #    end

        #end

        #@testset "all sites are covered, only once" begin
        #    for N = 2:10, nranks = 1:min(N, 5), parity = 0:1
        #        covered = zeros(Int, N, N, N)

        #        stu_ranges = get_all_ranges_stu(N, nranks, parity)
        #        @testset for it = 1:N, iu = 1:it, is = 1:N
        #            if (it + iu + is) % 2 == parity
        #                for (isrange, itrange, iurange) in stu_ranges
        #                    if (is in isrange && it in itrange && iu in iurange)
        #                        covered[is, it, iu] += 1
        #                    end
        #                end
        #                @test covered[is, it, iu] == 1
        #            end
        #        end
        #    end
        #end
    end
end
