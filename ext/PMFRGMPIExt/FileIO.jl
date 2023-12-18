using OrdinaryDiffEq, DiffEqCallbacks, PMFRG

function PMFRG.setCheckpoint(
    Directory::String,
    State,
    saved_values,
    Lam,
    Par,
    checkPointList,
    ParallelizationScheme::PMFRG.UseMPI,
)
    if ParallelizationScheme.rank == 0
        PMFRG.setCheckpoint(
            Directory,
            State,
            saved_values,
            Lam,
            Par,
            checkPointList,
            PMFRG.MultiThreaded(),
        )
    end
end

function PMFRG.saveMainOutput(
    Filename::String,
    Solution::ODESolution,
    saved_values::DiffEqCallbacks.SavedValues,
    Par::PMFRG.PMFRGParams,
    Group::String,
    ParallelizationScheme::PMFRG.UseMPI,
)
    if ParallelizationScheme.rank == 0
        PMFRG.saveMainOutput(
            Filename,
            Solution,
            saved_values,
            Par,
            Group,
            PMFRG.MultiThreaded(),
        )
    end
end

function PMFRG.saveCurrentState(
    DirPath::String,
    State::AbstractArray,
    saved_Values::DiffEqCallbacks.SavedValues,
    Lam::Real,
    Par::PMFRG.PMFRGParams,
    ParallelizationScheme::UseMPI,
)
    if ParallelizationScheme.rank == 0
        PMFRG.saveCurrentState(
            DirPath,
            State,
            saved_Values,
            Lam,
            Par,
            PMFRG.MultiThreaded(),
        )
    end
end

PMFRG.saveCurrentState(
    ::Nothing,
    ::AbstractArray,
    ::DiffEqCallbacks.SavedValues,
    ::Real,
    ::PMFRG.PMFRGParams,
    ::UseMPI,
) = nothing




function PMFRG.SetCompletionCheckmark(DirPath::String, ParallelizationScheme::PMFRG.UseMPI)
    if ParallelizationScheme.rank == 0
        PMFRG.SetCompletionCheckmark(DirPath, PMFRG.MultiThreaded())
    end
end

PMFRG.SetCompletionCheckmark(::Nothing, ::PMFRG.UseMPI) = nothing

function PMFRG.setupDirectory(
    DirPath,
    Par,
    ParallelizationScheme::PMFRG.UseMPI;
    overwrite = false,
)
    if ParallelizationScheme.rank == 0
        PMFRG.setupDirectory(DirPath, Par, PMFRG.MultiThreaded(); overwrite = overwrite)
    end
end
