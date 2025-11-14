using PMFRGCore, MPI, TimerOutputs, PencilArrays
include("mpi/MPI_Detail.jl")
import .MPI_Detail

function PMFRGCore.rebuildStateStruct!(
    StateStruct::PMFRGCore.StateType,
    PartialStateVector::PencilArray,
    setup,
)

    (; StateMPIBuff) = setup
    gatherState!(StateMPIBuff, PartialStateVector)
    PMFRGCore.unpackStateVector!(StateStruct, StateMPIBuff.data)
end



function PMFRGCore.getXBubble!(Workspace, Lam, ::PMFRGCore.UseMPI)
    (; X, Par) = Workspace
    (; N, np_vec) = Par.NumericalParams

    if MPI.Initialized()
        nranks = MPI.Comm_size(MPI.COMM_WORLD)
        rank = MPI.Comm_rank(MPI.COMM_WORLD)
        Lam_root = MPI.bcast(Lam, 0, MPI.COMM_WORLD)
        @assert (Lam_root == Lam) begin
            "Lambda differs between MPI ranks! " *
            "on rank 0: $Lam_root," *
            "on rank $rank: $Lam"
        end


        # if (ns+nt+nu)%2 == 0	# skip unphysical bosonic frequency combinations
        #
        if (3 * np_vec[1]) % 2 == 0
            # then 1,1,1 , with parity 1, is not right
            parity = 1
        else
            parity = 0
        end

        @timeit_debug "get_ranges" all_ranges = MPI_Detail.get_all_ranges_st(N, nranks)
        iurange_full = 1:N
        isrange, itranges = all_ranges[rank+1]
        itrange = vcat(itranges...)
        @timeit_debug "partition" PMFRG.getXBubblePartition!(
            X,
            Workspace.State,
            Workspace.Deriv,
            Par,
            Workspace.Buffer,
            Lam,
            isrange,
            itrange,
            iurange_full,
        )


        @timeit_debug "communication" for root = 0:(nranks-1)
            isrange_root, itranges_root = all_ranges[root+1]

            for itrange_root in itranges_root
                iurange_restrict = 1:itrange_root.stop
                iurange_abc = Par.Options.usesymmetry ? iurange_restrict : iurange_full

                MPI.Bcast!(
                    (@view X.a[:, isrange_root, itrange_root, iurange_abc]),
                    root,
                    MPI.COMM_WORLD,
                )
                MPI.Bcast!(
                    (@view X.b[:, isrange_root, itrange_root, iurange_abc]),
                    root,
                    MPI.COMM_WORLD,
                )
                MPI.Bcast!(
                    (@view X.c[:, isrange_root, itrange_root, iurange_abc]),
                    root,
                    MPI.COMM_WORLD,
                )

                MPI.Bcast!(
                    (@view X.Ta[:, isrange_root, itrange_root, iurange_full]),
                    root,
                    MPI.COMM_WORLD,
                )
                MPI.Bcast!(
                    (@view X.Tb[:, isrange_root, itrange_root, iurange_full]),
                    root,
                    MPI.COMM_WORLD,
                )
                MPI.Bcast!(
                    (@view X.Tc[:, isrange_root, itrange_root, iurange_full]),
                    root,
                    MPI.COMM_WORLD,
                )
                MPI.Bcast!(
                    (@view X.Td[:, isrange_root, itrange_root, iurange_full]),
                    root,
                    MPI.COMM_WORLD,
                )
            end
        end
    else
        @warn "MPI package used but not initialized" maxlog = 1
        PMFRGCore.getXBubblePartition!(X, State, Deriv, Par, Buffer, Lam, 1:N, 1:N, 1:N)
    end
end

function PMFRGCore.getChi(LocalState::PencilArray, Lam::Real, Par::PMFRGCore.PMFRGParams)
    GlobalState = gatherState(LocalState, Par)
    PMFRGCore.getChi(GlobalState, Lam, Par)
end
