"""
This module provides partitioning logic
for the MPI partitioning of the computation
of the arrays a,b,c (and Ta,Tb,Tc,Td)
in PMFRG.getXBubblePartition!

The ranges in the s,t indices are computed.

The S - T domain is partitioned using iteration pairing for load balancing.
Each rank works on frequencies near both ends of the spectrum (low and high),
or on a contiguous middle section, to balance the non-uniform computational costs
across the frequency range.
"""
module MPI_Detail
include("./partition.jl")
include("./best_partition_pairs.jl")
using .BestPartitionPairs: get_all_ranges_st

export get_all_ranges_st
end
