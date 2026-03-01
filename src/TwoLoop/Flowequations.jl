function getDerivTwoLoop!(Deriv, State, setup, Lam)
    X, Y, Buffs, Par = setup #use pre-allocated X and XTilde to reduce garbage collector time
    Workspace = TwoLoopWorkspace(Deriv, State, X, Y, Buffs, Par)

    getDFint!(Workspace, Lam)
    get_Self_Energy!(Workspace, Lam)
    getXBubble!(Workspace, Lam)

    symmetrizeBubble!(Workspace.X, Par)

    addToVertexFromBubble!(Workspace.Deriv.Gamma, Workspace.X)
    getYBubble!(Workspace, Lam)

    symmetrizeBubble!(Workspace.Y, Par)
    addToVertexFromBubble!(Workspace.Deriv.Gamma, Workspace.Y)

    symmetrizeVertex!(Workspace.Deriv.Gamma, Par)
    return
end

function getYBubble!(Workspace::TwoLoopWorkspace, Lam)
    (; X, Y, State, Par, Buffer) = Workspace
    getProp! = constructPropagatorFunction(Workspace, Lam)
    (; State, Par) = Workspace
    setZero!(Y)
    addToLeft2PartBubble!(Y, X, X, State.Gamma, getProp!, Par, Buffer)
    addToRight2PartBubble!(Y, X, X, State.Gamma, getProp!, Par, Buffer)
end

@inline function convertFreqArgsXT(ns, nt, nu, Nw)
    # @assert (ns+nt+nu) %2 != 0 "trying to convert wrong freqs $ns + $nt +  $nu = $(ns+nt+nu)"
    swapsites = ns * nu < 0
    ns, nt, nu = abs.((ns, nt, nu))
    ns = min(ns, Nw - 1 - (ns + Nw - 1) % 2)
    nt = min(nt, Nw - 1 - (nt + Nw - 1) % 2)
    nu = min(nu, Nw - 1 - (nu + Nw - 1) % 2)
    return ns, nt, nu, swapsites
end

"""
Computes a two-particle bubble in the s-Channel given two four-point functions (i.e. vertices or other bubbles).
Gamma is assumed to be a vertex.
To allow the computation using just the left part of a bubble, specification of the transpose is needed i.e. Transpose(X_L) = X_R and Transpose(X) = X, where X = XL + XR is the full bubble.
"""
function addTo2PartBubble!(
    ResultBubble::BubbleType,
    X,
    XTransp,
    Gamma::VertexType,
    getProp!::Function,
    addBTilde_Func!::Function,
    addB_Func!::Function,
    Par,
    Buffer,
)
    (; N, lenIntw, np_vec) = Par.NumericalParams
    @sync begin
        for is = 1:N, it = 1:N
            Threads.@spawn begin

                PropsXY = take!(Buffer.PropsXY)# get pre-allocated thread-safe buffers
                VBuffer = take!(Buffer.Vertex)
                XBuffer = take!(Buffer.X)
                ns = np_vec[is]
                nt = np_vec[it]
                for nw = (-lenIntw):(lenIntw-1) # Matsubara sum
                    spropX, spropY = getProp!(PropsXY, nw, nw + ns)
                    for iu = 1:N
                        nu = np_vec[iu]
                        if (ns + nt + nu) % 2 == 0# skip unphysical bosonic frequency combinations
                            continue
                        end
                        addBTilde_Func!(
                            ResultBubble,
                            X,
                            XTransp,
                            Gamma,
                            is,
                            it,
                            iu,
                            nw,
                            Par,
                            spropY,
                        )
                        if (!Par.Options.usesymmetry || nu <= nt)
                            addB_Func!(
                                ResultBubble,
                                X,
                                XTransp,
                                Gamma,
                                is,
                                it,
                                iu,
                                nw,
                                Par,
                                spropX,
                                VBuffer,
                                XBuffer,
                            )
                        end
                    end
                end
                put!(Buffer.PropsXY, PropsXY)
                put!(Buffer.Vertex, VBuffer)
                put!(Buffer.X, XBuffer)
            end
        end
    end
    return ResultBubble
end

@inline function compute2PartBubble!(
    ResultBubble::BubbleType,
    X,
    XTransp,
    Gamma::VertexType,
    getProp!::Function,
    addBTilde_Func!::Function,
    addB_Func!::Function,
    Par,
    Buffer,
)
    setZero!(ResultBubble)
    addTo2PartBubble!(
        ResultBubble,
        X,
        XTransp,
        Gamma,
        getProp!,
        addBTilde_Func!,
        addB_Func!,
        Par,
        Buffer,
    )
    symmetrizeBubble!(ResultBubble, Par)
end

@inline addToLeft2PartBubble!(
    ResultBubble::BubbleType,
    X,
    XTransp,
    Gamma::VertexType,
    getProp!::Function,
    Par,
    Buffer,
) = addTo2PartBubble!(
    ResultBubble,
    X,
    XTransp,
    Gamma,
    getProp!,
    addBLTilde!,
    addBL!,
    Par,
    Buffer,
)

@inline addToRight2PartBubble!(
    ResultBubble::BubbleType,
    X,
    XTransp,
    Gamma::VertexType,
    getProp!::Function,
    Par,
    Buffer,
) = addTo2PartBubble!(
    ResultBubble,
    X,
    XTransp,
    Gamma,
    getProp!,
    addBRTilde!,
    addBR!,
    Par,
    Buffer,
)

@inline computeLeft2PartBubble!(
    ResultBubble::BubbleType,
    X,
    XTransp,
    Gamma::VertexType,
    getProp!::Function,
    Par,
    Buffer,
) = compute2PartBubble!(
    ResultBubble,
    X,
    XTransp,
    Gamma,
    getProp!,
    addBLTilde!,
    addBL!,
    Par,
    Buffer,
)

@inline computeRight2PartBubble!(
    ResultBubble::BubbleType,
    X,
    XTransp,
    Gamma::VertexType,
    getProp!::Function,
    Par,
    Buffer,
) = compute2PartBubble!(
    ResultBubble,
    X,
    XTransp,
    Gamma,
    getProp!,
    addBRTilde!,
    addBR!,
    Par,
    Buffer,
)
