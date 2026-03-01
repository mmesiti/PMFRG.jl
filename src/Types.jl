"""
Abstract Param struct to dispatch between different PMFRG methods (i.e. Two-Loop or Parquet)
Assumed to have at least the fields:
System::Geometry
NumericalParams::NumericalParams
Options
"""
abstract type PMFRGParams end
abstract type PMFRGWorkspace end
abstract type AbstractParallelizationScheme end

struct MultiThreaded <: AbstractParallelizationScheme end
struct UseMPI <: AbstractParallelizationScheme end

"""
Struct to hold all relevant quantities that are needed throughout the computation. 

    N::Int = 24
    NiSigma::Int = 100 #Number of iSigma frequencies
    VDims::NTuple{4,Int} = (Npairs,N,N,N)
    accuracy::Float64 = 1e-6
    Lam_min::Float64 = exp(-10.),
    Lam_max::Float64 = exp(10.),
    usesymmetry::Bool = true
    MinimalOutput::Bool = false
    ex_freq::Float64 = (2*N-1)*pi*T
    lenIntw::Int = N
    lenIntw_acc::Int = 2*maximum((N,NiSigma,lenIntw)) # more accurate for less demanding sums
    np_vec::Array{Int,1} = collect(0:N-1)
    np_vec_iSigma::Array{Int,1} = collect(0:NiSigma-1)

"""
struct NumericalParams{F<:AbstractFloat}
    T::F # Temperature
    N::Int # Number of positive frequencies
    NiSigma::Int  # Number of positive iSigma frequencies

    accuracy::F # 
    Lam_min::F
    Lam_max::F

    lenIntw::Int
    lenIntw_acc::Int  # more accurate for less demanding sums
    np_vec::Array{Int,1}
    np_vec_iSigma::Array{Int,1}

    ex_freq::F
end

"""Converts other arguments to type of first argument. If that is not a float, converts everything to Float64."""
_convertToFloat(Primary::T, args...) where {T<:AbstractFloat} =
    convert.(T, (Primary, args...))
_convertToFloat(Primary::Real, args...) = convert.(Float64, (Primary, args...))

function NumericalParams(;
    T::Real = 0.5, # Temperature
    N::Integer = 24,
    NiSigma::Integer = N, #Number of iSigma frequencies
    accuracy::AbstractFloat = 1e-6, # convert type to float type
    Lam_min::AbstractFloat = exp(-10.0),
    Lam_max::AbstractFloat = exp(10.0),
    lenIntw::Int = N,
    lenIntw_acc::Int = 2 * maximum((N, NiSigma, lenIntw)), # more accurate for less demanding sums
    np_vec::Array{Int,1} = collect(0:(N-1)),
    np_vec_iSigma::Array{Int,1} = collect(0:(NiSigma-1)),
    ex_freq = (2 * N - 1) * pi * T,
    kwargs...,
)

    T, Lam_min, Lam_max, accuracy, ex_freq =
        _convertToFloat(T, Lam_min, Lam_max, accuracy, ex_freq)

    return NumericalParams(
        T,
        N,
        NiSigma,
        accuracy,
        Lam_min,
        Lam_max,
        lenIntw,
        lenIntw_acc,
        np_vec,
        np_vec_iSigma,
        ex_freq,
    )
end
abstract type AbstractOptions end
getPMFRGMethod(::Val{1}) = OneLoop()
getPMFRGMethod(::Val{2}) = TwoLoop()
getPMFRGMethod(::Val{n}) where {n} = MultiLoop(n)
getPMFRGMethod(n::Int) = getPMFRGMethod(Val(n))

abstract type AbstractVertex{T} end

abstract type AbstractSigma{T} end

struct HeisenbergSigma{T} <: AbstractSigma{T}
    data::Array{T,2}
end

HeisenbergSigma(NUnique::Int, NiSigma::Int, type = Float64) =
    HeisenbergSigma(zeros(type, NUnique, NiSigma))

Base.size(s::HeisenbergSigma) = size(s.data)
Base.size(s::HeisenbergSigma, i::Int) = size(s.data, i)
Base.IndexStyle(::Type{<:HeisenbergSigma}) = IndexLinear()
Base.getindex(s::HeisenbergSigma, i...) = getindex(s.data, i...)
Base.setindex!(s::HeisenbergSigma, val, i...) = setindex!(s.data, val, i...)
Base.fill!(s::HeisenbergSigma, x) = (fill!(s.data, x); s)
Base.copy(s::HeisenbergSigma) = HeisenbergSigma(copy(s.data))
Base.axes(s::HeisenbergSigma, i::Int) = axes(s.data, i)
Base.length(s::HeisenbergSigma) = length(s.data)
Base.eltype(::HeisenbergSigma{T}) where {T} = T
Base.eachindex(s::HeisenbergSigma) = eachindex(s.data)
Base.eachindex(a::HeisenbergSigma, b::HeisenbergSigma) = eachindex(a.data, b.data)
Base.similar(s::HeisenbergSigma, ::Type{T}, dims::Dims) where {T} =
    HeisenbergSigma(similar(s.data, T, dims))
Base.similar(s::HeisenbergSigma) = HeisenbergSigma(similar(s.data))
Base.maximum(s::HeisenbergSigma) = maximum(s.data)
Base.maximum(f, s::HeisenbergSigma; kwargs...) = maximum(f, s.data; kwargs...)

"""Struct containing information about four-point vertices"""
struct HeisenbergVertex{T} <: AbstractVertex{T}
    a::Array{T,4}
    b::Array{T,4}
    c::Array{T,4}
end


function HeisenbergVertex(VDims::Tuple, type)
    return HeisenbergVertex(
        zeros(type, VDims), # Gamma_a
        zeros(type, VDims), # Gamma_b
        zeros(type, VDims), # Gamma_c
    )
end
getVDims(Par::PMFRGParams) =
    (Par.System.Npairs, Par.NumericalParams.N, Par.NumericalParams.N, Par.NumericalParams.N)
HeisenbergVertex(Par::PMFRGParams) = HeisenbergVertex(getVDims(Par), _getFloatType(Par))

function setToBareVertex!(Gammac::AbstractArray{T,4}, couplings::AbstractVector) where {T}
    for Rj in eachindex(couplings, axes(Gammac, 1))
        Gammac[Rj, :, :, :] .= -couplings[Rj]
    end
    return Gammac
end

function setToBareVertex!(Gamma::HeisenbergVertex, couplings::AbstractVector)
    Gamma.a .= 0.0
    Gamma.b .= 0.0
    setToBareVertex!(Gamma.c, couplings)
    return Gamma
end

setToBareVertex!(Gamma::HeisenbergVertex, Par) =
    setToBareVertex!(Gamma, Par.System.couplings)

"""Allocates a frequency-dependent vertex type and initializes it with the bare vertex"""
BareVertex_Freq(Par) = setToBareVertex!(HeisenbergVertex(Par), Par.System.couplings)

"""Stores all information about the bare vertex (without frequencies)"""
struct BareHeisenbergVertex{T}
    c::Vector{T}
end

BareHeisenbergVertex(Par::PMFRGParams) =
    BareHeisenbergVertex(-convert.(_getFloatType(Par), Par.System.couplings))

"""Struct containing information about the (physical) ODE State, i.e. vertices"""
struct StateType{T}
    f_int::Array{T,1}
    iSigma::HeisenbergSigma{T}
    Gamma::HeisenbergVertex{T}
end

function StateType(NUnique::Int, NiSigma::Int, VDims::Tuple, type = Float64::Type)
    return StateType(
        zeros(type, NUnique), # fint
        HeisenbergSigma(NUnique, NiSigma, type), # iSigma
        HeisenbergVertex(VDims, type),
    )
end

StateType(Par::PMFRGParams) = StateType(
    Par.System.NUnique,
    Par.NumericalParams.NiSigma,
    getVDims(Par),
    _getFloatType(Par),
)

StateType(f_int, iSigma::Array{T,2}, Gamma::HeisenbergVertex) where {T} =
    StateType(f_int, HeisenbergSigma(iSigma), Gamma)
StateType(f_int, iSigma::HeisenbergSigma, Gammaa, Gammab, Gammac) =
    StateType(f_int, iSigma, HeisenbergVertex(Gammaa, Gammab, Gammac))
StateType(f_int, iSigma::Array{T,2}, Gammaa, Gammab, Gammac) where {T} =
    StateType(f_int, HeisenbergSigma(iSigma), HeisenbergVertex(Gammaa, Gammab, Gammac))


abstract type AbstractBubble{T} end

"""Struct storing information about Bubble functions, i.e. rhs derivatives of vertex flow equations."""
struct BubbleType{T} <: AbstractBubble{T}
    a::Array{T,4} #"a" type bubble
    b::Array{T,4}
    c::Array{T,4}

    Ta::Array{T,4} #"a-Tilde" type bubble
    Tb::Array{T,4}
    Tc::Array{T,4}
    Td::Array{T,4}
end

function BubbleType(VDims::Tuple, type = Float64)
    return BubbleType(
        zeros(type, VDims),
        zeros(type, VDims),
        zeros(type, VDims),
        zeros(type, VDims),
        zeros(type, VDims),
        zeros(type, VDims),
        zeros(type, VDims),
    )
end

const VertexOrBubble = Union{StateType,BubbleType,HeisenbergVertex}

BubbleType(Par::PMFRGParams) = BubbleType(getVDims(Par), _getFloatType(Par))

function constructBubbleFromVertex!(B::BubbleType, Gamma::HeisenbergVertex)
    B.a .= Gamma.a
    B.b .= Gamma.b
    B.c .= Gamma.c
    for Rij in axes(Gamma.a, 1),
        s in axes(Gamma.a, 2),
        t in axes(Gamma.a, 3),
        u in axes(Gamma.a, 4)

        B.Ta[Rij, s, t, u] = -Gamma.a[Rij, t, s, u]
        B.Tb[Rij, s, t, u] = -Gamma.c[Rij, t, s, u]
        B.Tc[Rij, s, t, u] = -Gamma.b[Rij, t, s, u]
        B.Td[Rij, s, t, u] = Gamma.c[Rij, t, u, s]
    end
    return B
end
constructBubbleFromVertex(Gamma::HeisenbergVertex) =
    constructBubbleFromVertex!(BubbleType(size(Gamma.a)), Gamma)

"""Struct containing the observables that are saved at every step"""
struct Observables{T}
    Chi::Vector{T}
    iSigma::Matrix{T}
    f_int::Vector{T}
    MaxVa::Vector{T}
    MaxVb::Vector{T}
    MaxVc::Vector{T}
end

struct VertexBufferType{T}
    Va12::Vector{T}
    Vb12::Vector{T}
    Vc12::Vector{T}

    Va34::Vector{T}
    Vb34::Vector{T}
    Vc34::Vector{T}

    Vc21::Vector{T}
    Vc43::Vector{T}
end
VertexBufferType(type, Npairs) = VertexBufferType((zeros(type, Npairs) for _ = 1:8)...)
##
RecursiveArrayTools.ArrayPartition(x::StateType) =
    ArrayPartition(x.f_int, x.iSigma.data, x.Gamma.a, x.Gamma.b, x.Gamma.c)
StateType(Arr::ArrayPartition) = StateType(Arr.x...)
