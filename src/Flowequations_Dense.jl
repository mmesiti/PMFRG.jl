function getDeriv!(Deriv,State,setup::Tuple{BubbleType,T,OneLoopParams},Lam) where T
    (X,Buffs,Par) = setup #use pre-allocated X and XTilde to reduce garbage collector time
    Workspace = OneLoopWorkspace(Deriv,State,X,Buffs,Par)

    getDFint!(Workspace,Lam)
    get_Self_Energy!(Workspace,Lam)

    getXBubble!(Workspace,Lam)

    symmetrizeBubble!(Workspace.X,Par)

    addToVertexFromBubble!(Workspace.Deriv.Γ,Workspace.X)
    symmetrizeVertex!(Workspace.Deriv.Γ,Par)
    flush(stdout)
    return
end

function getDerivVerbose!(Deriv,State,setup,Lam)
    (X,Buffs,Par) = setup #use pre-allocated X and XTilde to reduce garbage collector time
    print("Workspace:\n\t") 
    @time Workspace = OneLoopWorkspace(Deriv,State,X,Buffs,Par)
    print("getDFint:\n\t") 
    @time getDFint!(Workspace,Lam)
    print("get_Self_Energy:\n\t") 
    @time get_Self_Energy!(Workspace,Lam)
    print("getVertexDeriv:\n\t") 
    @time getXBubble!(Workspace,Lam)
    print("Symmetry:\n\t") 
    @time begin
        symmetrizeBubble!(Workspace.X,Par)
        addToVertexFromBubble!(Workspace.Deriv.Γ,Workspace.X)
        symmetrizeVertex!(Workspace.Deriv.Γ,Par)
    end
    flush(stdout)
    return
end


function getDFint!(Workspace::PMFRGWorkspace,Lam::double)
    @unpack State,Deriv,Par = Workspace 
    @unpack T,lenIntw_acc = Par.NumericalParams 
    NUnique = Par.System.NUnique 
	
	@inline γ(x,nw) = gamma_(State.γ,x,nw)
	@inline iG(x,nw) = iG_(State.γ,x,Lam,nw,T)
	@inline iS(x,nw) = iS_(State.γ,x,Lam,nw,T)

	Theta(Lam,w) = w^2/(w^2+Lam^2)
	
	for x in 1:NUnique
		sumres = 0.
		for nw in -lenIntw_acc:lenIntw_acc-1
			w = get_w(nw,T)
			sumres += iS(x,nw)/iG(x,nw)*Theta(Lam,w) *γ(x,nw)/w
		end
		Deriv.f_int[x] = -3/2*T*sumres
	end
end

"""
Computes a single-particle (i.e. self-energy) bubble. Allows specification of function type, i.e. what vertices are used since this is different if a bubble function is inserted as opposed to a vertex.
"""
function addTo1PartBubble!(Dgamma::AbstractArray,XT1_::Function,XT2_::Function,Prop,Par)

    @unpack T,N,Ngamma,lenIntw_acc,np_vec_gamma = Par.NumericalParams
    @unpack siteSum,invpairs,Nsum,OnsitePairs = Par.System

	Threads.@threads for iw1 in 1:Ngamma
		nw1 = np_vec_gamma[iw1]
    	for (x,Rx) in enumerate(OnsitePairs)
			for nw in -lenIntw_acc: lenIntw_acc-1
				jsum = 0.
				w1pw = nw1+nw+1 #w1 + w: Adding two Matsubara frequencies gives a +1
				w1mw = nw1-nw
				for k_spl in 1:Nsum[Rx]
					@unpack m,ki,xk = siteSum[k_spl,Rx]
					jsum += (XT1_(ki,w1pw,0,w1mw)+2*XT2_(ki,w1pw,0,w1mw))*Prop(xk,nw)*m
				end
				Dgamma[x,iw1] += -T *jsum #For the self-energy derivative, the factor of 1/2 must be in the propagator
			end
		end
	end
    return Dgamma
end

"""
Computes a single-particle (i.e. self-energy) bubble. Can only be used if B is a bubble function
"""
function addTo1PartBubble!(Dgamma::AbstractArray,X::BubbleType,Prop,Par)
	invpairs = Par.System.invpairs
	@inline XTa_(Rij,s,t,u) = XT_(X.Ta,X.Ta,Rij,s,t,u,invpairs[Rij],Par.NumericalParams.N)
	@inline XTc_(Rij,s,t,u) = XT_(X.Tc,X.Tc,Rij,s,t,u,invpairs[Rij],Par.NumericalParams.N)
    addTo1PartBubble!(Dgamma,XTa_,XTc_,Prop,Par)
end

"""
Computes a single-particle (i.e. self-energy) bubble. Can only be used if argument is a vertex
"""
function addTo1PartBubble!(Dgamma::AbstractArray,Γ::VertexType,Prop,Par)
    invpairs = Par.System.invpairs
    # @warn "addTo1PartBubble! for vertices is not tested yet!"
	@inline Γa_(Rij,s,t,u) = V_(Γ.a,invpairs[Rij],t,u,s,Rij,Par.NumericalParams.N) # Tilde-type can be obtained by permutation of vertices
	@inline Γb_(Rij,s,t,u) = V_(Γ.b,invpairs[Rij],t,u,s,Rij,Par.NumericalParams.N) # cTilde corresponds to b type vertex!
    addTo1PartBubble!(Dgamma,Γa_,Γb_,Prop,Par)
end

function compute1PartBubble!(Dgamma::AbstractArray,ΓorX,Prop,Par)
	setZero!(Dgamma)
    addTo1PartBubble!(Dgamma,ΓorX,Prop,Par)
end
function get_Self_Energy!(Workspace::PMFRGWorkspace,Lam) 
	Par = Workspace.Par
	@inline iS(x,nw) = iS_(Workspace.State.γ,x,Lam,nw,Par.NumericalParams.T)/2
	compute1PartBubble!(Workspace.Deriv.γ,Workspace.State.Γ,iS,Par)
end
# @inline getXBubble!(Workspace::PMFRGWorkspace,Lam) = getXBubble!(Workspace,Lam,Val(Workspace.Par.System.NUnique)) 

function getXBubble!(Workspace::PMFRGWorkspace,Lam)
	Par = Workspace.Par
    @unpack T,N,lenIntw,np_vec = Par.NumericalParams 
    PropsBuffers = Workspace.Buffer.Props 
    VertexBuffers = Workspace.Buffer.Vertex
	 
	iG(x,nw) = iG_(Workspace.State.γ,x,Lam,nw,T)
	iSKat(x,nw) = iSKat_(Workspace.State.γ,Workspace.Deriv.γ,x,Lam,nw,T)

	function getKataninProp!(BubbleProp,nw1,nw2)
		# BubbleProp = SpinFRGLattices.MMatrix{NUnique,NUnique,Float64,NUnique*NUnique}(undef)
		for i in 1:Par.System.NUnique, j in 1:Par.System.NUnique
			BubbleProp[i,j] = iSKat(i,nw1) *iG(j,nw2)* T
		end
		return SpinFRGLattices.SMatrix(BubbleProp)
		# SpinFRGLattices.MMatrix[iSKat(i,nw1) *iG(j,nw2)* T for i in 1:NUnique, j in 1:NUnique]
	end
	@sync begin
		for is in 1:N,it in 1:N
			Threads.@spawn begin
				BubbleProp = PropsBuffers[Threads.threadid()] # get pre-allocated thread-safe buffers
				Buffer = VertexBuffers[Threads.threadid()]
				ns = np_vec[is]
				nt = np_vec[it]
				for iu in 1:N
					nu = np_vec[iu]
					if (ns+nt+nu)%2 == 0	# skip unphysical bosonic frequency combinations
						continue
					end
					for nw in -lenIntw:lenIntw-1 # Matsubara sum
						sprop = getKataninProp!(BubbleProp,nw,nw+ns) 
						addXTilde!(Workspace,is,it,iu,nw,sprop) # add to XTilde-type bubble functions
						if(!Par.Options.usesymmetry || nu<=nt)
							addX!(Workspace,is,it,iu,nw,sprop,Buffer)# add to X-type bubble functions
						end
					end
				end
			end
		end
	end
end

@inline function mixedFrequencies(ns,nt,nu,nwpr)
	nw1=Int((ns+nt+nu-1)/2)
    nw2=Int((ns-nt-nu-1)/2)
    nw3=Int((-ns+nt-nu-1)/2)
    nw4=Int((-ns-nt+nu-1)/2)
	wpw1 = nwpr + nw1+1
    wpw2 = nwpr + nw2+1
    wmw3 = nwpr - nw3
    wmw4 = nwpr - nw4
	# @assert (ns + wmw3 +wmw4)%2 != 0 "error in freq"
	return wpw1,wpw2,wmw3,wmw4
end

"""
adds part of X functions in Matsubara sum at nwpr containing the site summation for a set of s t and u frequencies. This is the most numerically demanding part!
"""
function addX!(Workspace::PMFRGWorkspace, is::Integer, it::Integer, iu::Integer, nwpr::Integer, Props,Buffer)
	@unpack State,X,Par = Workspace 
	@unpack Va12,Vb12,Vc12,Va34,Vb34,Vc34,Vc21,Vc43 = Buffer 
	@unpack N,np_vec = Par.NumericalParams
	@unpack Npairs,Nsum,siteSum,invpairs = Par.System
	Npairs = Workspace.Par.System.Npairs
	ns = np_vec[is]
	nt = np_vec[it]
	nu = np_vec[iu]
	wpw1,wpw2,wmw3,wmw4 = mixedFrequencies(ns,nt,nu,nwpr)

	bufferV_!(Va12, State.Γ.a , ns, wpw1, wpw2, invpairs,N)
	bufferV_!(Vb12, State.Γ.b , ns, wpw1, wpw2, invpairs,N)
	bufferV_!(Vc12, State.Γ.c , ns, wpw1, wpw2, invpairs,N)

	bufferV_!(Va34, State.Γ.a , ns, wmw3, wmw4, invpairs,N)
	bufferV_!(Vb34, State.Γ.b , ns, wmw3, wmw4, invpairs,N)
	bufferV_!(Vc34, State.Γ.c , ns, wmw3, wmw4, invpairs,N)
	
	bufferV_!(Vc21, State.Γ.c , ns, wpw2, wpw1, invpairs,N)
	bufferV_!(Vc43, State.Γ.c , ns, wmw4, wmw3, invpairs,N)
	# get fields of siteSum struct as Matrices for better use of LoopVectorization
	S_ki = siteSum.ki
	S_kj = siteSum.kj
	S_xk = siteSum.xk
	S_m = siteSum.m

	for Rij in 1:Npairs
		#loop over all left hand side inequivalent pairs Rij
		Xa_sum = 0. #Perform summation on this temp variable before writing to State array as Base.setindex! proved to be a bottleneck!
		Xb_sum = 0.
		Xc_sum = 0.
		@turbo unroll = 1 for k_spl in 1:Nsum[Rij]
			#loop over all Nsum summation elements defined in geometry. This inner loop is responsible for most of the computational effort! 
			# @unpack ki,kj,m,xk = siteSum[k_spl,Rij]
			ki,kj,m,xk = S_ki[k_spl,Rij],S_kj[k_spl,Rij],S_m[k_spl,Rij],S_xk[k_spl,Rij]
			Ptm = Props[xk,xk]*m

			Xa_sum += (
				+Va12[ki] * Va34[kj] 
				+Vb12[ki] * Vb34[kj] * 2
			)* Ptm

			Xb_sum += (
				+Va12[ki] * Vb34[kj]
				+Vb12[ki] * Va34[kj]
				+Vb12[ki] * Vb34[kj]
			)* Ptm
			
			Xc_sum += (
				+Vc12[ki] * Vc34[kj]
				+Vc21[ki] * Vc43[kj]
			)* Ptm
		end
		X.a[Rij,is,it,iu] += Xa_sum
		X.b[Rij,is,it,iu] += Xb_sum
		X.c[Rij,is,it,iu] += Xc_sum
    end
    return
end
##
function addXTilde!(Workspace::PMFRGWorkspace, is::Integer, it::Integer, iu::Integer, nwpr::Integer, Props)

	@unpack State,X,Par = Workspace 
	@unpack N,np_vec = Par.NumericalParams
	@unpack Npairs,invpairs,PairTypes,OnsitePairs = Par.System

	@inline Va_(Rij,s,t,u) = V_(State.Γ.a,Rij,s,t,u,invpairs[Rij],N)
	@inline Vb_(Rij,s,t,u) = V_(State.Γ.b,Rij,s,t,u,invpairs[Rij],N)
	@inline Vc_(Rij,s,t,u) = V_(State.Γ.c,Rij,s,t,u,invpairs[Rij],N)
	ns = np_vec[is]
	nt = np_vec[it]
	nu = np_vec[iu]
	wpw1,wpw2,wmw3,wmw4 = mixedFrequencies(ns,nt,nu,nwpr)

	#Xtilde only defined for nonlocal pairs Rij != Rii
	for Rij in 1:Npairs
		Rij in OnsitePairs  && continue
		#loop over all left hand side inequivalent pairs Rij
		Rji = invpairs[Rij] # store pair corresponding to Rji (easiest case: Rji = Rij)
		@unpack xi,xj = PairTypes[Rij]

		#These values are used several times so they are saved locally
		Va12 = Va_(Rji, wpw1, ns, wpw2)
		Va21 = Va_(Rij, wpw2, ns, wpw1)
		Va34 = Va_(Rji, wmw3, ns, wmw4)
		Va43 = Va_(Rij, wmw4, ns, wmw3)

		Vb12 = Vb_(Rji, wpw1, ns, wpw2)
		Vb21 = Vb_(Rij, wpw2, ns, wpw1)
		Vb34 = Vb_(Rji, wmw3, ns, wmw4)
		Vb43 = Vb_(Rij, wmw4, ns, wmw3)

		Vc12 = Vc_(Rji, wpw1, ns, wpw2)
		Vc21 = Vc_(Rij, wpw2, ns, wpw1)
		Vc34 = Vc_(Rji, wmw3, ns, wmw4)
		Vc43 = Vc_(Rij, wmw4, ns, wmw3)

	    X.Ta[Rij,is,it,iu] += (
			(+Va21 * Va43
			+2*Vc21 * Vc43) * Props[xi,xj]
			+(Va12 * Va34
			+2*Vc12 * Vc34)* Props[xj,xi]
		)
		
	    X.Tb[Rij,is,it,iu] += (
			(+Va21 * Vc43
			+Vc21 * Vc43
			+Vc21 * Va43) * Props[xi,xj]

			+(Va12 * Vc34
			+Vc12 * Vc34
			+Vc12 * Va34)* Props[xj,xi]
		)
		Vb12 = Vb_(Rji, wpw1, wpw2, ns)
		Vb21 = Vb_(Rij, wpw2, wpw1, ns)
		Vb34 = Vb_(Rji, wmw3, wmw4, ns)
		Vb43 = Vb_(Rij, wmw4, wmw3, ns)

		Vc12 = Vc_(Rji, wpw1, wpw2, ns)
		Vc21 = Vc_(Rij, wpw2, wpw1, ns)
		Vc34 = Vc_(Rji, wmw3, wmw4, ns)
		Vc43 = Vc_(Rij, wmw4, wmw3, ns)


	    X.Tc[Rij,is,it,iu] += (
			(+Vb21 * Vb43
			+Vc21 * Vc43
			) * Props[xi,xj]
			+(Vb12 * Vb34
			+Vc12 * Vc34
	    	)* Props[xj,xi]
		)
    end
end
const SingleElementMatrix = Union{SpinFRGLattices.SMatrix{1,1},SpinFRGLattices.MMatrix{1,1}}

"""Use multiple dispatch to treat the common special case in which the propagator does not depend on site indices to increase performance"""
function addXTilde!(Workspace::PMFRGWorkspace, is::Integer, it::Integer, iu::Integer, nwpr::Integer, Props::SingleElementMatrix)

	@unpack State,X,Par = Workspace 
	@unpack N,np_vec = Par.NumericalParams
	@unpack Npairs,invpairs,PairTypes,OnsitePairs = Par.System

	@inline Va_(Rij,s,t,u) = V_(State.Γ.a,Rij,s,t,u,invpairs[Rij],N)
	@inline Vb_(Rij,s,t,u) = V_(State.Γ.b,Rij,s,t,u,invpairs[Rij],N)
	@inline Vc_(Rij,s,t,u) = V_(State.Γ.c,Rij,s,t,u,invpairs[Rij],N)
	ns = np_vec[is]
	nt = np_vec[it]
	nu = np_vec[iu]
	wpw1,wpw2,wmw3,wmw4 = mixedFrequencies(ns,nt,nu,nwpr)
	Prop = only(Props)
	#Xtilde only defined for nonlocal pairs Rij != Rii
	for Rij in 1:Npairs
		Rij in OnsitePairs  && continue
		#loop over all left hand side inequivalent pairs Rij
		Rji = invpairs[Rij] # store pair corresponding to Rji (easiest case: Rji = Rij)
		@unpack xi,xj = PairTypes[Rij]

		#These values are used several times so they are saved locally
		Va12 = Va_(Rji, wpw1, ns, wpw2)
		Va21 = Va_(Rij, wpw2, ns, wpw1)
		Va34 = Va_(Rji, wmw3, ns, wmw4)
		Va43 = Va_(Rij, wmw4, ns, wmw3)

		Vb12 = Vb_(Rji, wpw1, ns, wpw2)
		Vb21 = Vb_(Rij, wpw2, ns, wpw1)
		Vb34 = Vb_(Rji, wmw3, ns, wmw4)
		Vb43 = Vb_(Rij, wmw4, ns, wmw3)

		Vc12 = Vc_(Rji, wpw1, ns, wpw2)
		Vc21 = Vc_(Rij, wpw2, ns, wpw1)
		Vc34 = Vc_(Rji, wmw3, ns, wmw4)
		Vc43 = Vc_(Rij, wmw4, ns, wmw3)

	    X.Ta[Rij,is,it,iu] += (
			+Va21 * Va43
			+2*Vc21 * Vc43
			+Va12 * Va34
			+2*Vc12 * Vc34
		)*Prop
		
	    X.Tb[Rij,is,it,iu] += (
			+Va21 * Vc43
			+Vc21 * Vc43
			+Vc21 * Va43

			+Va12 * Vc34
			+Vc12 * Vc34
			+Vc12 * Va34
		)*Prop
		Vb12 = Vb_(Rji, wpw1, wpw2, ns)
		Vb21 = Vb_(Rij, wpw2, wpw1, ns)
		Vb34 = Vb_(Rji, wmw3, wmw4, ns)
		Vb43 = Vb_(Rij, wmw4, wmw3, ns)

		Vc12 = Vc_(Rji, wpw1, wpw2, ns)
		Vc21 = Vc_(Rij, wpw2, wpw1, ns)
		Vc34 = Vc_(Rji, wmw3, wmw4, ns)
		Vc43 = Vc_(Rij, wmw4, wmw3, ns)


	    X.Tc[Rij,is,it,iu] += (
			+Vb21 * Vb43
			+Vc21 * Vc43
			+Vb12 * Vb34
			+Vc12 * Vc34
		)*Prop
    end
end


"""Use multiple dispatch to treat the common special case in which the propagator does not depend on site indices to increase performance"""
function addX!(Workspace::PMFRGWorkspace, is::Integer, it::Integer, iu::Integer, nwpr::Integer, Props::SingleElementMatrix,Buffer)
	@unpack State,X,Par = Workspace 
	@unpack Va12,Vb12,Vc12,Va34,Vb34,Vc34,Vc21,Vc43 = Buffer 
	@unpack N,np_vec = Par.NumericalParams
	@unpack Npairs,Nsum,siteSum,invpairs = Par.System
	Npairs = Workspace.Par.System.Npairs

	ns = np_vec[is]
	nt = np_vec[it]
	nu = np_vec[iu]
	wpw1,wpw2,wmw3,wmw4 = mixedFrequencies(ns,nt,nu,nwpr)

	bufferV_!(Va12, State.Γ.a , ns, wpw1, wpw2, invpairs,N)
	bufferV_!(Vb12, State.Γ.b , ns, wpw1, wpw2, invpairs,N)
	bufferV_!(Vc12, State.Γ.c , ns, wpw1, wpw2, invpairs,N)

	bufferV_!(Va34, State.Γ.a , ns, wmw3, wmw4, invpairs,N)
	bufferV_!(Vb34, State.Γ.b , ns, wmw3, wmw4, invpairs,N)
	bufferV_!(Vc34, State.Γ.c , ns, wmw3, wmw4, invpairs,N)
	
	bufferV_!(Vc21, State.Γ.c , ns, wpw2, wpw1, invpairs,N)
	bufferV_!(Vc43, State.Γ.c , ns, wmw4, wmw3, invpairs,N)
	# get fields of siteSum struct as Matrices for better use of LoopVectorization
	S_ki = siteSum.ki
	S_kj = siteSum.kj
	S_m = siteSum.m
	Prop = only(Props)
	# Prop = Props
	for Rij in 1:Npairs
		#loop over all left hand side inequivalent pairs Rij
		Xa_sum = 0. #Perform summation on this temp variable before writing to State array as Base.setindex! proved to be a bottleneck!
		Xb_sum = 0.
		Xc_sum = 0.
		@turbo unroll = 1 for k_spl in 1:Nsum[Rij]
			#loop over all Nsum summation elements defined in geometry. This inner loop is responsible for most of the computational effort! 
			# @unpack ki,kj,m,xk = siteSum[k_spl,Rij]
			ki,kj,m = S_ki[k_spl,Rij],S_kj[k_spl,Rij],S_m[k_spl,Rij]

			Xa_sum += (
				+Va12[ki] * Va34[kj] 
				+Vb12[ki] * Vb34[kj] * 2
			)* m

			Xb_sum += (
				+Va12[ki] * Vb34[kj]
				+Vb12[ki] * Va34[kj]
				+Vb12[ki] * Vb34[kj]
			)* m
			
			Xc_sum += (
				+Vc12[ki] * Vc34[kj]
				+Vc21[ki] * Vc43[kj]
			)* m
		end
		X.a[Rij,is,it,iu] += Xa_sum*Prop
		X.b[Rij,is,it,iu] += Xb_sum*Prop
		X.c[Rij,is,it,iu] += Xc_sum*Prop
    end
	return
end

"""Use symmetries and identities to compute the rest of bubble functions"""
function symmetrizeBubble!(X::BubbleType,Par::PMFRGParams)
    N = Par.NumericalParams.N
    @unpack Npairs,OnsitePairs = Par.System
    usesymmetry = Par.Options.usesymmetry
    # use the u <--> t symmetry
    if(usesymmetry)
        Threads.@threads for it in 1:N
            for iu in it+1:N, is in 1:N, Rij in 1:Npairs
                X.a[Rij,is,it,iu] = -X.a[Rij,is,iu,it]
                X.b[Rij,is,it,iu] = -X.b[Rij,is,iu,it]
                X.c[Rij,is,it,iu] = (
                + X.a[Rij,is,it,iu]+
                - X.b[Rij,is,it,iu]+
                + X.c[Rij,is,iu,it])
            end
        end
    end
    #local definitions of X.Tilde vertices
    for iu in 1:N, it in 1:N, is in 1:N, R in OnsitePairs
        X.Ta[R,is,it,iu] = X.a[R,is,it,iu]
        X.Tb[R,is,it,iu] = X.b[R,is,it,iu]
        X.Tc[R,is,it,iu] = X.c[R,is,it,iu]
        X.Td[R,is,it,iu] = -X.c[R,is,iu,it]
    end
    @. X.Td= X.Ta - X.Tb - X.Tc
end

function addToVertexFromBubble!(Γ::VertexType,X::BubbleType)
    Threads.@threads for iu in axes(Γ.a,4)
        for it in axes(Γ.a,3), is in axes(Γ.a,2), Rij in axes(Γ.a,1)
            Γ.a[Rij,is,it,iu] += X.a[Rij,is,it,iu] - X.Ta[Rij,it,is,iu] + X.Ta[Rij,iu,is,it]
            Γ.b[Rij,is,it,iu] += X.b[Rij,is,it,iu] - X.Tc[Rij,it,is,iu] + X.Tc[Rij,iu,is,it]
            Γ.c[Rij,is,it,iu] += X.c[Rij,is,it,iu] - X.Tb[Rij,it,is,iu] + X.Td[Rij,iu,is,it]
        end
    end
    return Γ
end


function symmetrizeVertex!(Γ::VertexType,Par)
	N = Par.NumericalParams.N
	for iu in 1:N, it in 1:N, is in 1:N, R in Par.System.OnsitePairs
		Γ.c[R,is,it,iu] = -Γ.b[R,it,is,iu]
	end
end

##
getChi(State::ArrayPartition, Lam::double,Par::PMFRGParams,Numax) = getChi(State.x[2],State.x[5], Lam,Par,Numax)
getChi(State::ArrayPartition, Lam::double,Par::PMFRGParams) = getChi(State.x[2],State.x[5], Lam,Par)

function getChi(gamma::AbstractArray,Γc::AbstractArray, Lam::double,Par::PMFRGParams,Numax)
	@unpack T,N,lenIntw_acc,np_vec = Par.NumericalParams
	@unpack Npairs,invpairs,PairTypes,OnsitePairs = Par.System

	@inline iG(x,w) = iG_(gamma,x, Lam,w,T)
	@inline Vc_(Rij,s,t,u) = V_(Γc,Rij,s,t,u,invpairs[Rij],N)

	Chi = zeros(Npairs,N)

	@inbounds Threads.@threads for Rij in 1:Npairs
		@unpack xi,xj = PairTypes[Rij]
		for i_nu in 1:Numax
       		n_nu = np_vec[i_nu]
		
			for nK in -lenIntw_acc:lenIntw_acc-1
				if Rij in OnsitePairs
					Chi[Rij,i_nu] += T * iG(xi,nK) * iG(xi,nK+n_nu)
				end
				for nK2 in -lenIntw_acc:lenIntw_acc-1
					npwpw2 = n_nu+nK+nK2+1
					wmw2 = nK-nK2
					#use that Vc_0 is calculated from Vb
					GGGG = iG(xi,nK)*iG(xi,nK+n_nu) * iG(xj,nK2)*iG(xj,nK2+n_nu)
					Chi[Rij,i_nu] += T^2 * GGGG *Vc_(Rij,n_nu,npwpw2,wmw2)
                end
            end
        end
    end
	return(Chi)
end

function getChi(gamma::AbstractArray,Γc::AbstractArray, Lam::Real,Par::PMFRGParams)
	@unpack T,N,lenIntw_acc,np_vec = Par.NumericalParams
	@unpack Npairs,invpairs,PairTypes,OnsitePairs = Par.System

	@inline iG(x,w) = iG_(gamma,x, Lam,w,T)
	@inline Vc_(Rij,s,t,u) = V_(Γc,Rij,s,t,u,invpairs[Rij],N)

	Chi = zeros(Npairs)

	@inbounds Threads.@threads for Rij in 1:Npairs
		@unpack xi,xj = PairTypes[Rij]
		for nK in -lenIntw_acc:lenIntw_acc-1
			if Rij in OnsitePairs
				Chi[Rij,1] += T * iG(xi,nK) ^2
			end
			for nK2 in -lenIntw_acc:lenIntw_acc-1
				npwpw2 = nK+nK2+1
				wmw2 = nK-nK2
				#use that Vc_0 is calculated from Vb
				GGGG = iG(xi,nK)^2 * iG(xj,nK2)^2
				Chi[Rij] += T^2 * GGGG *Vc_(Rij,0,npwpw2,wmw2)
			end
        end
    end
	return(Chi)
end