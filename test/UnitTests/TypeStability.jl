
"""gets empty Workspace."""
function getEmptyWorkspace(Par::PMFRG.OneLoopParams)
    State = PMFRG.InitializeState(Par)
    Deriv = similar(State)
    setup = PMFRG.AllocateSetup(Par)
    (X, Buffs, Par) = setup
    PMFRG.OneLoopWorkspace(Deriv, State, X, Buffs, Par)
end

convertToSMatrix(M::PMFRG.MMatrix{N,N,T,NN}) where {N,T,NN} = PMFRG.SMatrix{N,N,T,NN}(M)
convertToSVector(V::PMFRG.MVector{N,T}) where {N,T} = PMFRG.SVector{N,T}(V)

"""Checks whether vertex derivative has type instabilities"""
function test_OneLoopAllocations(Par = Params(getPolymer(2)))
    WS = getEmptyWorkspace(Par)
    (PropsXBuff, PropsYBuff) = take!(WS.Buffer.PropsXY)
    VBuff = take!(WS.Buffer.Vertex)

    # Static array versions
    SPropsX = convertToSVector(PropsXBuff)
    SPropsY = convertToSMatrix(PropsYBuff)
    test_OneLoopAllocations(Par, WS, SPropsX, SPropsY, VBuff)

    # Plain array versions
    PropsYMatrix = Matrix(SPropsY)
    test_OneLoopAllocations(Par, WS, SPropsX, PropsYMatrix, VBuff)

    put!(WS.Buffer.PropsXY, (PropsXBuff, PropsYBuff))
    put!(WS.Buffer.Vertex, VBuff)
end
function test_OneLoopAllocations(Par, WS, PropsX, PropsY, VBuff)
    PMFRG.addX!(WS.X, WS.State, WS.Par, 1, 1, 2, 2, PropsX, VBuff) # compile functions
    PMFRG.addXTilde!(WS.X, WS.State, WS.Par, 1, 1, 2, 2, PropsY) # compile functions

    aX = @allocated PMFRG.addX!(WS.X, WS.State, WS.Par, 1, 1, 2, 2, PropsX, VBuff)
    aXT = @allocated PMFRG.addXTilde!(WS.X, WS.State, WS.Par, 1, 1, 2, 2, PropsY)
    test_TypeStability(Par, WS, aX, aXT)
end

function test_TypeStability(Par, WS, aX, aXT)
    @testset verbose = true "Allocations and Type stability" begin
        @testset "Type concreteness of Parameters " begin
            @test isconcretetype(typeof(Par))
        end
        @testset "Type concreteness of Workspace " begin
            @test isconcretetype(typeof(WS))
        end
        @testset "zero allocations of addX! " begin
            @test aX == 0
        end
        @testset "zero allocations of add XTilde! " begin
            @test aXT == 0
        end
    end
end

function getEmptyWorkspace(Par::PMFRG.TwoLoopParams)
    State = PMFRG.InitializeState(Par)
    Deriv = similar(State)
    setup = PMFRG.AllocateSetup(Par)
    (X, Y, Buffs, Par) = setup
    PMFRG.TwoLoopWorkspace(Deriv, State, X, Y, Buffs, Par)
end
function test_TwoLoopAllocations(Par = Params(getPolymer(2), TwoLoop()))
    WS = getEmptyWorkspace(Par)
    PropB_0X, PropB_0Y = take!(WS.Buffer.PropsXY)

    VB = take!(WS.Buffer.Vertex)
    XB = take!(WS.Buffer.X)

    PropBX = convertToSVector(PropB_0X)
    PropBY = convertToSMatrix(PropB_0Y)

    PropBXVector = Vector(PropBX)
    PropBYMatrix = Matrix(PropBY)

    test_TwoLoopAllocations(Par, WS, PropBX, PropBY, VB, XB)
    test_TwoLoopAllocations(Par, WS, PropBXVector, PropBYMatrix, VB, XB)

    put!(WS.Buffer.Vertex, VB)
    put!(WS.Buffer.X, XB)
    put!(WS.Buffer.PropsXY, (PropB_0X, PropB_0Y))

end
function test_TwoLoopAllocations(Par, WS, PropBX, PropBY, VB, XB)
    PMFRG.addBL!(WS.Y, WS.X, WS.X, WS.State.Gamma, 1, 1, 2, 2, Par, PropBX, VB, XB) #compile functions
    PMFRG.addBR!(WS.Y, WS.X, WS.X, WS.State.Gamma, 1, 1, 2, 2, Par, PropBX, VB, XB) #compile functions
    PMFRG.addBLTilde!(WS.Y, WS.X, WS.X, WS.State.Gamma, 1, 1, 2, 2, Par, PropBY) #compile functions
    PMFRG.addBRTilde!(WS.Y, WS.X, WS.X, WS.State.Gamma, 1, 1, 2, 2, Par, PropBY) #compile functions

    aX = @allocated PMFRG.addBL!(
        WS.Y,
        WS.X,
        WS.X,
        WS.State.Gamma,
        1,
        1,
        2,
        2,
        Par,
        PropBX,
        VB,
        XB,
    )
    aX += @allocated PMFRG.addBR!(
        WS.Y,
        WS.X,
        WS.X,
        WS.State.Gamma,
        1,
        1,
        2,
        2,
        Par,
        PropBX,
        VB,
        XB,
    )
    aXT = @allocated PMFRG.addBLTilde!(
        WS.Y,
        WS.X,
        WS.X,
        WS.State.Gamma,
        1,
        1,
        2,
        2,
        Par,
        PropBY,
    )
    aXT += @allocated PMFRG.addBRTilde!(
        WS.Y,
        WS.X,
        WS.X,
        WS.State.Gamma,
        1,
        1,
        2,
        2,
        Par,
        PropBY,
    )

    test_TypeStability(Par, WS, aX, aXT)

end
# Todo: try to write unit tests using @code_warntype  

function capture_stdout(f)
    stdout_orig = stdout
    (rd, wr) = redirect_stdout()
    f()
    close(wr)
    redirect_stdout(stdout_orig)
    read(rd, String)
end
