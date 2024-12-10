# from deprecated PolynomialMatrices.jl, cudos to Neveritt

using Polynomials, Compat, DSP, FFTW, DataStructures

import Polynomials: coeffs, degree
import Base: size
import LinearAlgebra: det

const SymbolLike = Union{Symbol,AbstractString,Char}
const ForwardOrdering = Base.Order.ForwardOrdering

struct PolyMatrix{T,M,W,N} <: AbstractArray{Polynomials.Polynomial{T},N}
    coeffs::SortedDict{Int,M,ForwardOrdering}
    dims::NTuple{N,Int}
    
    @compat function (::Type{PolyMatrix})(
        coeffs::SortedDict{Int,M,ForwardOrdering}, dims::NTuple{N,Int}, ::Type{Val{W}}) where {M,N,W}
      T = eltype(M)
      _truncate!(coeffs, dims, T)
      new{T,M,Val{W},N}(coeffs, dims)
    end
  end
  
  function truncate!(p::PolyMatrix{T,M,Val{V},N},
    ϵ=Base.rtoldefault(T)*length(p)*degree(p)) where {T,M,V,N}
    _truncate!(coeffs(p), size(p), T, ϵ)
  end
  
  function _truncate!(coeffs::SortedDict{Int,M,ForwardOrdering},
    dims::NTuple{N,Int}, ::Type{T},
    ϵ::Real=Base.rtoldefault(real(T))*prod(dims)*length(coeffs)) where {T,M,N}

    v1 = coeffs[first(keys(coeffs))]
    for (st,k,v) in semitokens(coeffs)
      nonz = true
      for i in eachindex(v)
        v[i]  = abs(v[i]) < ϵ ? zero(v[i]) : v[i]
        nonz &= abs(v[i]) < ϵ
      end
      if nonz
        delete!((coeffs,st)) # all entries are zero => remove semitoken
      end
    end
    if length(coeffs) == zero(T)
      insert!(coeffs, 0, zero(v1))
    end
    coeffs
  end
    
  function PolyMatrix(PM::M1) where M1<:AbstractArray
    var = count(!iszero,PM) > 0 ? variable(PM[findfirst(x -> x != zero(x), PM)]) :
                         variable(Polynomial(T[]))   # default to Polynomials default variable
    PolyMatrix(PM, Val{@compat Symbol(var)})
  end
  
  function PolyMatrix(PM::AbstractArray{Polynomial{T},N}, ::Type{Val{W}}) where {T,N,W}
    N <= 2 || error("higher order arrays not supported at this point")
    M = typeof(similar(PM, T)) # NOTE: Is there a more memory-efficient way to obtain M?
    c = SortedDict(Dict{Int,M}())
    # find the union of all index sets of all polynomials in the matrix PM
    S = Set{Int}()
    for p in PM
      for elem in eachindex(coeffs(p))
        push!(S, elem)
      end
    end
    # initialize all elements to zero
    for idx in S
      insert!(c, idx-1, zeros(T, size(PM)...)) # TODO change to spzeros when Julia drops support for 0.4.7 (there are no sparse vectors)
    end
    # copy all elements
    for pidx in eachindex(PM)
      pc = coeffs(PM[pidx])
      for eidx in eachindex(pc)
        c[eidx-1][pidx] = pc[eidx]
      end
    end
    PolyMatrix(c, size(PM), Val{W})
  end
  

  size(p::PolyMatrix) = p.dims
  size(p::PolyMatrix, i::Int) = i ≤ length(p.dims) ? p.dims[i] : 1
  coeffs(p::PolyMatrix) = p.coeffs
  degree(p::PolyMatrix{T,M,Val{W},N}) where {T,M,W,N}  = last(coeffs(p))[1]

  function det(p::PolyMatrix{T,M,Val{W},N}) where {T,M,W,N}
    size(p,1) == size(p,2) || throw(DimensionMismatch("the polynomial matrix must be square for a determinant to be defined"))
    n  = size(p,1)
    dn = (degree(p))*n+1
    # copy all elements into three-dimensional matrix
    A = zeros(n,n,dn)
    for (k,v) in coeffs(p)
      A[:,:,k+1] = v
    end
    # take fft and evaluate determinant at each interpolation point
    B = fft(A,3)
    a = [det(B[:,:,k]) for k = 1:dn]
    # interpolate using fft
    ar = _truncate(T, ifft(a))
    return Polynomial(ar, W)
  end
  
  function _truncate(::Type{T}, a::AbstractArray{T2}) where {T<:Real,T2}
    r = similar(a, T)
    r = real(a)
  end
  
  function _truncate(::Type{T}, a::AbstractArray{T2}) where {T<:Integer,T2}
    r = similar(a, T)
    for i in eachindex(a)
      r[i] = convert(T, round(real(a[i])))
    end
    r
  end