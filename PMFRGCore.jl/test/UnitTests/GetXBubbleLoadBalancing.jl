include("../../src/GetXBubbleLoadBalancing.jl")


function getGetXBubblePartitionWorkloadInSingleIndexRange(
    isrange,
    itrange,
    iurange,
    singleindexrange,
)
    work = 0.0
    for i in singleindexrange
        is, it = getIsItFromSingleIndex(isrange, itrange, i)
        ns = is - 1
        nt = it - 1
        for iu in iurange
            nu = iu - 1
            if (ns + nt + nu) % 2 == 0# skip unphysical bosonic frequency combinations
                continue
            end
            work += addxtildework
            if (nu <= nt)
                work += addxwork
            end
        end
    end
    return work
end

"""
We should make sure we do not leave any index out
"""
function testGetXBubbleLoadBalancing()

    @testset "(is,it) -> single index conversion" begin
        @test getSingleIsItIndex(5:10, 5:10, 5, 5) == 1
        @test getSingleIsItIndex(5:10, 5:10, 10, 10) == 6 * 6

        @test getSingleIsItIndex(5:10, 5:10, 5, 6) == 2
        @test getSingleIsItIndex(5:10, 5:10, 6, 5) == 7
    end

    @testset "single index -> (is,it) conversion" begin
        @test getIsItFromSingleIndex(5:10, 5:10, 1) == (5, 5)
        @test getIsItFromSingleIndex(5:10, 5:10, 6 * 6) == (10, 10)
        @test getIsItFromSingleIndex(5:10, 5:10, 2) == (5, 6)
        @test getIsItFromSingleIndex(5:10, 5:10, 7) == (6, 5)

    end

    @testset "round trip" begin
        @testset for i = 1:36
            @test let (is, it) = getIsItFromSingleIndex(5:10, 5:10, i)
                getSingleIsItIndex(5:10, 5:10, is, it) == i
            end
        end
    end

    @testset "All indices covered, once" begin


        @testset for (isrange, itrange, iurange, nthreads) in
                     [(5:11, 6:11, 5:11, 5), (1:5, 1:5, 1:5, 32)]

            singleindexranges = getSingleindexranges(isrange, itrange, iurange, nthreads)
            @test length(singleindexranges) == nthreads

            @testset for i = 1:length(isrange)*length(itrange)
                @test sum(i in r for r in singleindexranges) == 1
            end
        end
    end
    @testset "load is balanced" begin
        isrange = 1:11
        itrange = 1:11
        iurange = 1:11

        nthreads = 5
        singleindexranges = getSingleindexranges(isrange, itrange, iurange, nthreads)

        each_thread_work = [
            getGetXBubblePartitionWorkloadInSingleIndexRange(
                isrange,
                itrange,
                iurange,
                singleindexrange,
            ) for singleindexrange in singleindexranges
        ]

        min_work = min(each_thread_work...)
        max_work = max(each_thread_work...)

        @test max_work - min_work < length(iurange) * (addxwork + addxtildework)
    end


end
