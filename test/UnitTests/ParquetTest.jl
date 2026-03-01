function test_SDE(
    Par::PMFRG.ParquetParams = Params(
        getPolymer(2),
        Parquet(),
        T = 1.0,
        N = 30,
        NiSigma = 30,
        accuracy = 1e-5,
        Lam_min = 0.0,
        Lam_max = 100.0,
        usesymmetry = false,
        lenIntw = 800,
        lenIntw_acc = 1000,
        maxIterSDE = 300,
        maxIterBSE = 80,
        MinimalOutput = true,
    );
    tol = 1e-5,
)
    Workspace = PMFRG.SetupParquet(Par)
    Lam = 0.0

    (; State, Gamma0, X, B0, BX, Par, Buffer) = Workspace
    @inline Prop(x, nw) = 1 / 6 * PMFRG.iG_(State.iSigma, x, Lam, nw, Par.NumericalParams.T)

    test_FirstBSEIteration!(Workspace, Lam)
    iSigma1 = copy(State.iSigma)
    iSigma2 = copy(State.iSigma)


    PMFRG.compute1PartBubble!(iSigma1, B0, Prop, Par)
    PMFRG.compute1PartBubble_BS!(iSigma2, State.Gamma, Gamma0, Prop, Par)
    @testset "SDE iSigma" begin
        for i in eachindex(iSigma1, iSigma2)
            @test iSigma1[i] ≈ iSigma2[i] atol = tol
        end
    end
end

function test_FirstBSEIteration!(Workspace, Lam)
    (; OldState, State, Gamma0, X, B0, BX, Par, Buffer) = Workspace

    getProp! = PMFRG.constructPropagatorFunction(Workspace, Lam)

    PMFRG.PMFRG.computeLeft2PartBubble!(
        B0,
        Gamma0,
        Gamma0,
        State.Gamma,
        getProp!,
        Par,
        Buffer,
    )
end

function test_SDE_FP(Par::PMFRG.ParquetParams)
    Workspace1 = PMFRG.SetupParquet(Par)
    Workspace2 = PMFRG.SetupParquet(Par)
    Lam = 0.0
    test_FirstBSEIteration!(Workspace1, Lam)
    test_FirstBSEIteration!(Workspace2, Lam)

    PMFRG.iterateSDE!(Workspace1, Lam)
    # PMFRG.iterateSDE_FP!(Workspace2,Lam)
    @testset "SDE FixedPoint" begin
        State1 = Workspace1.State
        State2 = Workspace2.State
        for i in eachindex(State1.iSigma, State2.iSigma)
            @test State1.iSigma[i] ≈ State2.iSigma[i]
        end
        for field in (:a, :b, :c)
            V1 = getproperty(State1.Gamma, field)
            V2 = getproperty(State2.Gamma, field)
            for i in eachindex(V1, V2)
                @test V1[i] ≈ V2[i]
            end
        end

    end
end


function test_iterate!(Workspace::PMFRG.ParquetWorkspace, Lam::Real; SDECheatfactor = 1)
    (; OldState, State, I, Gamma0, X, B0, BX, Par, Buffer) = Workspace

    maxIterBSE = Par.Options.maxIterBSE


    Tol_Vertex = 1E16 #Initial tolerance level

    getProp! = constructPropagatorFunction(Workspace, Lam)
    iter = 0


    PMFRG.writeTo!(OldState.Gamma, State.Gamma)

    # PMFRG.iterateSDE!(Workspace,Lam)

    PMFRG.computeLeft2PartBubble!(B0, Gamma0, Gamma0, State.Gamma, getProp!, Par, Buffer)
    PMFRG.computeLeft2PartBubble!(BX, X, X, State.Gamma, getProp!, Par, Buffer)

    PMFRG.getXFromBubbles!(X, B0, BX) #TODO when applicable, this needs to be generalized for beyond-Parquet approximations 
    PMFRG.getVertexFromChannels!(State.Gamma, I, X)
    PMFRG.symmetrizeVertex!(State.Gamma, Par)

    # PMFRG.iterateSDE!(Workspace,Lam,SDECheatfactor = SDECheatfactor)
    # PMFRG.iterateSDE_FP!(Workspace,Lam)
    PMFRG.iterateSDE_FP!(State.iSigma, State.Gamma, Gamma0, Lam, Par)
    Tol_Vertex = reldist(OldState.Gamma, State.Gamma)

    isnan(Tol_Vertex) && @warn ": BSE: Solution diverged after $iter iterations\n"
    println(maximum.((PMFRG.ArrayPartition(OldState), PMFRG.ArrayPartition(State))))
    return Workspace
end

function test_iterate_FP!(Workspace::PMFRG.ParquetWorkspace, Lam::Real)
    (; OldState, State, I, Gamma0, X, B0, BX, Par, Buffer) = Workspace

    maxIterBSE = Par.Options.maxIterBSE

    OldStateArr, StateArr = PMFRG.ArrayPartition.((OldState, State))

    function FixedPointFunction!(State_Arr, OldState_Arr)
        anyisnan(OldState_Arr) && return State_Arr
        iSigma = OldState_Arr.x[2]

        State = StateType(State_Arr.x...)
        OldState = StateType(OldState_Arr.x...)
        getProp! = constructPropagatorFunction(iSigma, Lam, Par)

        PMFRG.computeLeft2PartBubble!(
            B0,
            Gamma0,
            Gamma0,
            OldState.Gamma,
            getProp!,
            Par,
            Buffer,
        )
        PMFRG.computeLeft2PartBubble!(BX, X, X, OldState.Gamma, getProp!, Par, Buffer)

        PMFRG.getXFromBubbles!(X, B0, BX) #TODO when applicable, this needs to be generalized for beyond-Parquet approximations 
        PMFRG.getVertexFromChannels!(State.Gamma, I, X)
        PMFRG.symmetrizeVertex!(State.Gamma, Par)


        WS_State = PMFRG.ArrayPartition(Workspace.State)
        WS_State .= State_Arr
        # PMFRG.iterateSDE_FP!(Workspace,Lam)
        PMFRG.iterateSDE_FP!(State.iSigma, State.Gamma, Gamma0, Lam, Par)

        println(maximum.((OldState_Arr, State_Arr)))
        return State_Arr
    end

    # FixedPointFunction!(StateArr,OldStateArr)
    s = afps!(FixedPointFunction!, OldStateArr, iters = 1, vel = 0.0, ep = 1.0, tol = 0.0)
    # StateArr .= s.x
    # OldStateArr .= StateArr
    # isnan(Tol_Vertex) && @warn ": BSE: Solution diverged after $iter iterations\n"
    return Workspace
end

function test_BSE(Par::PMFRG.ParquetParams)
    Workspace1 = PMFRG.SetupParquet(Par)
    Workspace2 = PMFRG.SetupParquet(Par)
    Lam = 0.0

    for _ = 1:4
        test_iterate!(Workspace1, Lam)
        test_iterate_FP!(Workspace2, Lam)
        println("FP!:")
    end
    Gamma1 = Workspace1.State.Gamma
    Gamma2 = Workspace2.State.Gamma
    # @testset "BSE test" begin
    #     for k in fieldnames(HeisenbergVertex)
    #         G1,G2 = getfield.((Gamma1,Gamma2),Ref(k))
    #         for i in eachindex(G1,G2)
    #             @test G1[i] ≈ G2[i]
    #         end
    #     end
    # end
    println(maximum(Gamma1.c))
    println(maximum(Gamma2.c))
    println(maximum(Workspace1.State.iSigma))
    println(maximum(Workspace2.State.iSigma))
end
