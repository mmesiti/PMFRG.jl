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

    ranges = []
    start = 1
    total_work = 0
    for i = 1:length(isrange)*length(itrange)
        is, it = getIsItFromSingleIndex(isrange, itrange, i)
        ns = is - 1
        nt = it - 1
        optimal_work_this_thread = work_by_thread * (length(ranges) + 1)
        prediff = total_work - optimal_work_this_thread
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
        postdiff = total_work - optimal_work_this_thread
        if prediff <= 0 <= postdiff
            if postdiff > -prediff
                stop = i - 1
            else
                stop = i
            end
            push!(ranges, start:stop)
            start = stop + 1
        end
    end
    return ranges
end
