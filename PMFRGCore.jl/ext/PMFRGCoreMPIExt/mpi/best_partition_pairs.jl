module BestPartitionPairs
include("./partition.jl")
import .Partition: partitions

struct TooManyRanksError <: Exception
    rank_t::Int
    N::Int
end

function Base.show(io::IO, e::TooManyRanksError)
    println(io, "Too many ranks ($(e.rank_t)) for N=$(e.N)")
end



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

    if N < 2 * nranks_t
        throw(TooManyRanksError(nranks_t, N))
    end

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


function get_imbalance_from_ranges(ranges_per_all_ranks, N)
    loads_per_rank = Load[]

    for (isrange, itranges) in ranges_per_all_ranks
        rank_load = Load(0, 0)
        for is in isrange
            for itrange in itranges
                for it in itrange
                    # Simple model:
                    # - addXLoad increases linearly with it
                    # - addXTildeLoad is constant
                    # We neglect the effect of parity.

                    rank_load = rank_load + Load(it, N)
                end
            end
        end
        push!(loads_per_rank, rank_load)
    end

    if isempty(loads_per_rank)
        return 0.0
    end

    max_load_addX = maximum(l.addXLoad for l in loads_per_rank)
    mean_load_addX = sum(l.addXLoad for l in loads_per_rank) / length(loads_per_rank)
    imbalance_addX = (max_load_addX - mean_load_addX) / mean_load_addX

    max_load_addXTilde = maximum(l.addXTildeLoad for l in loads_per_rank)
    mean_load_addXTilde =
        sum(l.addXTildeLoad for l in loads_per_rank) / length(loads_per_rank)
    imbalance_addXTilde = (max_load_addXTilde - mean_load_addXTilde) / mean_load_addXTilde

    return imbalance_addX, imbalance_addXTilde
end

function _all_ns_x_nt_factorizations(nranks)
    [(nranks_s, div(nranks, nranks_s)) for nranks_s = 1:nranks if (nranks % nranks_s) == 0]
end

"""
Print load balancing factorization results in a readable format.
"""
function _print_load_balancing_results(
    N::Int,
    nranks::Int,
    all_ranges_best_X,
    all_ranges_best_XTilde,
    imbalance_X::Float64,
    imbalance_XTilde::Float64,
)
    println("="^70)
    println("MPI Load Balancing Factorization Results (N=$N, nranks=$nranks)")
    println("="^70)

    if all_ranges_best_X == all_ranges_best_XTilde
        # Both objectives agree on the same factorization
        nranks_s_best = length(unique(r[1] for r in all_ranges_best_X))
        nranks_t_best = div(nranks, nranks_s_best)
        println(
            "Optimal factorization: nranks_s × nranks_t = $nranks_s_best × $nranks_t_best",
        )
        println("  Relative imbalance (addX):      $(round(imbalance_X * 100, digits=2))%")
        println(
            "  Relative imbalance (addXTilde): $(round(imbalance_XTilde * 100, digits=2))%",
        )
    else
        # Different optimal factorizations for the two objectives
        println("Warning: Different optimal factorizations for addX and addXTilde")
        println()

        nranks_s_X = length(unique(r[1] for r in all_ranges_best_X))
        nranks_t_X = div(nranks, nranks_s_X)
        println("Best for addX: nranks_s × nranks_t = $nranks_s_X × $nranks_t_X")
        println("  Relative imbalance: $(round(imbalance_X * 100, digits=2))%")

        nranks_s_XTilde = length(unique(r[1] for r in all_ranges_best_XTilde))
        nranks_t_XTilde = div(nranks, nranks_s_XTilde)
        println(
            "Best for addXTilde: nranks_s × nranks_t = $nranks_s_XTilde × $nranks_t_XTilde",
        )
        println("  Relative imbalance: $(round(imbalance_XTilde * 100, digits=2))%")
        println()
    end
    println("="^70)
end

"""
For given values of N and nranks, returns the ranges in s and t
for all ranks using the pairing-based load balancing approach.

- Splits along the s direction in equal parts
  differing at most by unity
- Splits only along the t direction using iteration pairing

Returns: Vector of Tuples, where each tuple contains:
  (isrange, (itrange1, [itrange2]))
The it component is a tuple of 1 or 2 ranges depending on whether
the ranges touch.
"""
function get_all_ranges_st(N::Int, nranks::Int)
    imbalance_X = Inf
    imbalance_XTilde = Inf
    all_ranges_best_XTilde =
        Vector{Tuple{UnitRange{Int64},Tuple{Vararg{UnitRange{Int64}}}}}()
    all_ranges_best_X = Vector{Tuple{UnitRange{Int64},Tuple{Vararg{UnitRange{Int64}}}}}()


    for (nranks_s, nranks_t) in _all_ns_x_nt_factorizations(nranks)

        try
            all_ranges = Vector{Tuple{UnitRange{Int64},Tuple{Vararg{UnitRange{Int64}}}}}()

            s_partitions = partitions(N, nranks_s)
            for rank_s = 0:(nranks_s-1)

                isrange = s_partitions[1+rank_s]

                for rank_t = 0:(nranks_t-1)
                    itranges = _get_ranges_t(N, nranks_t, rank_t)
                    push!(all_ranges, (isrange, itranges))
                end
            end

            candidate_imbalance_X, candidate_imbalance_XTilde =
                get_imbalance_from_ranges(all_ranges, N)
            if candidate_imbalance_X < imbalance_X
                all_ranges_best_X = all_ranges
                imbalance_X = candidate_imbalance_X
            end
            if candidate_imbalance_XTilde < imbalance_XTilde
                all_ranges_best_XTilde = all_ranges
                imbalance_XTilde = candidate_imbalance_XTilde
            end

        catch e
            if e isa TooManyRanksError
                # Skip this factorization as it's invalid (N < 2*nranks_t)
                continue
            else
                rethrow(e)
            end
        end

    end


    _print_load_balancing_results(
        N,
        nranks,
        all_ranges_best_X,
        all_ranges_best_XTilde,
        imbalance_X,
        imbalance_XTilde,
    )

    println("Using addX-optimal factorization")
    all_ranges_best_X
end

end
