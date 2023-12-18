using RecursiveArrayTools, PMFRG

function PMFRG.writeOutput(
    State::ArrayPartition,
    saved_values,
    Lam,
    Par,
    ParallelizationScheme::PMFRG.UseMPI,
)
    if ParallelizationScheme.rank == 0
        PMFRG.writeOutput(State, saved_values, Lam, Par, PMFRG.MultiThreaded())
    end
end
