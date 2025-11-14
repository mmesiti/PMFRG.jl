module BestPartitionPairs

struct TooManyRanksError <: Exception
    rank_tu::Int
    N::Int
end

function Base.show(io::IO, e::TooManyRanksError)
    println(io, "Too many ranks ($(e.rank_tu)) for N=$(e.N)")
end

"""Returns the t and u ranges for a given rank out of nranks,
   trying to produce the most balanced partition
   (where the weight is the number of sites of the given parity
   in the triangle where iu <= it).
"""
function _get_ranges_t(N::Int, nranks_t::Int, rank_t::Int)
    fences(N, npieces) = [round(Int, i * N / npieces, RoundDown) for i = 0:npieces]
    starts(N, npieces) = 1 .+ fences(N, npieces)[1:(end-1)]
    ends(N, npieces) = fences(N, npieces)[2:end]
    partitions(N, npieces) = [s:e for (s, e) in zip(starts(N, npieces), ends(N, npieces))]

    all_ranges = partitions(N, nranks_t * 2)

    low_range = all_ranges[rank_t+1]
    high_range = all_ranges[nranks_t*2-rank_t]
    tlow_start, tlow_stop = low_range.start, low_range.stop
    thigh_start, thigh_stop = high_range.start, high_range.stop


    if tlow_stop < thigh_start - 1
        (tlow_start:tlow_stop, thigh_start:thigh_stop)
    else
        (tlow_start:thigh_stop,)
    end

end


function get_imbalance_from_ranges(N::Int, nranks::Int, ranges_per_all_ranks, parity::Int)
    min_nsites = typemax(Int64)
    max_nsites = 0
    for (irank, (isrange, itrange, e)) in enumerate(ranges_per_all_ranks)
        nsites = _count_sites(itrange, iurange, isrange, parity)
        min_nsites = (nsites < min_nsites) ? nsites : min_nsites
        max_nsites = (nsites > max_nsites) ? nsites : max_nsites
    end
    return (max_nsites - min_nsites) / max_nsites
end



"""Returns an estimate of the load imbalance among the ranks,
   in the range 0.0-1.0.
 """
function get_imbalance(N, nranks, get_ranges_func, parity)
    all_ranges = get_ranges_func(N, nranks, parity)
    get_imbalance_from_ranges(N, nranks, all_ranges, parity)
end


function _all_ns_x_ntu_factorizations(nranks)
    [(nranks_s, div(nranks, nranks_s)) for nranks_s = 1:nranks if (nranks % nranks_s) == 0]
end

"""
For given values of N, nranks and parity choices,
returns the ranges in s,t,u that correspond to the best balance partition,
for all the ranks.
"""
function get_all_ranges_stu(N::Int, nranks::Int, parity::Int)
    imbalance = 1.0

    all_ranges_best = Vector{Tuple{UnitRange{Int64},UnitRange{Int64},UnitRange{Int64}}}()

    for (nranks_s, nranks_tu) in _all_ns_x_ntu_factorizations(nranks)
        all_ranges = Vector{Tuple{UnitRange,UnitRange,UnitRange}}()
        for rank = 0:(nranks-1)

            rank_s = rank % nranks_s
            range_s = partitions(N, nranks_s)[1+rank_s]

            rank_tu = div(rank, nranks_s, RoundToZero)
            range_tu = _get_ranges_tu(N, nranks_tu, rank_tu)

            range_stu = (range_s, range_tu...)
            push!(all_ranges, range_stu)
        end
        candidate_imbalance = get_imbalance_from_ranges(N, nranks, all_ranges, parity)
        if candidate_imbalance < imbalance
            all_ranges_best = all_ranges
            imbalance = candidate_imbalance
        end
    end
    all_ranges_best
end

"""Returns the number of sites having iu<=it<=Ntu"""
function get_number_of_sites_in_triangle(Ntu)
    Int(Ntu * (Ntu + 1) / 2)
end

"""Returns the number of sites having iu<=it<=Ntu with a given parity."""
function _get_number_of_sites_eo(Ntu)
    total_elements = get_number_of_sites_in_triangle(Ntu)
    even = if (Ntu % 2 == 0)
        nhalf = Ntu / 2
        2 * nhalf * (nhalf + 1) / 2
    else
        nhalf = (Ntu - 1) / 2
        2 * nhalf * (nhalf + 1) / 2 + nhalf + 1
    end

    odd = total_elements - even
    even, odd
end
end
