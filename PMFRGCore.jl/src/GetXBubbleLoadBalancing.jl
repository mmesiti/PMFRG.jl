"""
Auxiliary functions to find best work split
between different threads
for the current implementation of getXBubblePartition!.

This means, partitioning the is-range and it-range in a way
so that the work is equal.
"""

"""First of all we need a function that unrolls
the is- and it- range into a single index."""

function getSingleIsItIndex(isrange, itrange, is, it)
    soffset = isrange.start
    toffset = itrange.start

    tstride = length(itrange)

    return (it - toffset + 1) + tstride * (is - soffset)
end

"""
The reverse
"""
function getIsItFromSingleIndex(isrange, itrange, singleindex)
    @assert 0 < singleindex <= length(isrange) * length(itrange)
    soffset = isrange.start
    toffset = itrange.start

    tstride = length(itrange)

    it = (singleindex - 1) % tstride + toffset
    is = Int(floor((singleindex - 1) / tstride)) + soffset
    return is, it

end


function getAllWork(isrange, itrange, iurange)
    total_work = 0.0
    addxtildework = 0.1
    addxwork = 1.0

    for i = 1:length(isrange)*length(itrange)
        is, it = getIsItFromSingleIndex(isrange, itrange, i)
        ns = is - 1
        nt = it - 1
        for iu in iurange
            nu = iu - 1
            if (ns + nt + nu) % 2 == 0# skip unphysical bosonic frequency combinations
                continue
            end
            total_work += addxtildework
            if (nu <= nt)
                total_work += addxwork
            end
        end
    end
    return total_work
end

addxtildework = 0.1
addxwork = 1.0


function getSingleindexranges(isrange, itrange, iurange, nthreads)
    allwork = getAllWork(isrange, itrange, iurange)
    work_by_thread = allwork / nthreads

    function work_by_i(i, isrange, itrange)
        @assert 0 < i <= length(isrange) * length(itrange)
        is, it = getIsItFromSingleIndex(isrange, itrange, i)
        w = 0
        for iu in iurange
            if (is + it + iu - 1) % 2 == 0# skip unphysical bosonic frequency combinations
                continue
            end
            w += addxtildework
            if (iu <= it)
                w += addxwork
            end
        end
        return w
    end
    work_by_threadid_cumulative(threadid) = threadid * work_by_thread

    # find matches
    range_stops = [0]
    i = 0
    cumulative_work_by_i = 0
    for threadid = 1:nthreads
        pre_diff = cumulative_work_by_i - work_by_threadid_cumulative(threadid)
        while cumulative_work_by_i < work_by_threadid_cumulative(threadid)
            if cumulative_work_by_i ≈ work_by_threadid_cumulative(threadid)
                break
            else
                i += 1
                cumulative_work_by_i += work_by_i(i, isrange, itrange)
            end
        end
        post_diff = cumulative_work_by_i - work_by_threadid_cumulative(threadid)

        if abs(pre_diff) < abs(post_diff)
            cumulative_work_by_i -= work_by_i(i, isrange, itrange)
            i -= 1
        end
        push!(range_stops, i)
    end
    ranges = [(range_stops[i]+1):range_stops[i+1] for i = 1:nthreads]

    return ranges
end
