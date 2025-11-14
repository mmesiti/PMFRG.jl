module BestPartitionPairs

"""Returns a pair of t ranges for a given rank out of nranks,
   trying to assign the same number of t-slices to all ranks.
   The pair will be composed of two ranges
   of the same length when possible
   (or of length differing by 1 at most)
   whose distance from the average value of t
   is equal (or it could differ by 1).
   
   When the two ranges touch, only one range is returned
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

struct Load
    addXLoad::Int64
    addXTildeLoad::Int64
end
import Base: +
+(a::Load, b::Load) = Load(a.addXLoad + b.addXLoad, a.addXTildeLoad + b.addXTildeLoad)



function get_imbalance_from_ranges(N::Int, nranks::Int, ranges_per_all_ranks, parity::Int)
    min_nsites = typemax(Int64)
    max_nsites = 0
    for (irank, (isrange, itranges)) in enumerate(ranges_per_all_ranks)
        nsites = _count_sites(itranges, isrange, parity)
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

"""
For given values of N and nranks, returns the ranges in s and t
for all ranks using the pairing-based load balancing approach.

Simplified version:
- Always returns 1:N for isrange (no splitting in s direction)
- Splits only along the t direction using iteration pairing
- Parity is not considered 

Returns: Vector of Tuples, where each tuple contains:
  (isrange, (itrange1, [itrange2]))
The it component is a tuple of 1 or 2 ranges depending on whether
the ranges touch.
"""
function get_all_ranges_st(N::Int, nranks::Int)
    all_ranges = Vector{Tuple{UnitRange{Int64},Tuple{Vararg{UnitRange{Int64}}}}}()

    for rank = 0:(nranks-1)
        isrange = 1:N  # Full range for s, no splitting yet
        itranges = _get_ranges_t(N, nranks, rank)
        push!(all_ranges, (isrange, itranges))
    end

    return all_ranges
end

end
