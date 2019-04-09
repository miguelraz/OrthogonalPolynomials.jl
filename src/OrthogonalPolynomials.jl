module OrthogonalPolynomials

# 0. Our mission
# Dear codes, should you choose to accept it, is to learn and implement the following
# mathematical in Julia, whilst learning the magic of multiple dispatch as a design
# paradigm. Ready? Steady? Go!
# If you are a beginner and this code just looks like squiggles, worry not!
# There is a video tutorial online at my youtube channel, BrainRPG.
# Happy Hacking!
# - Miguel Raz Guzman Macedo
# http://people.math.sfu.ca/~cbm/aands/page_789.jpg

# 1. First design - Laguerre only
# Note the _names with underscore are used as convention for package internals.
_d(n, α=0) = binomial(n+α, n)
_b(n, m) = n-m+1
_c(α, m) = m*(α+m)
_k(n, α=0) = [-_b(n,i)*inv(_c(α, i)) for i in 1:n]
_f(x) = x
_α = 0

_a(x, n, ks = _k(n,α), i=0) = i == n ? :(1) : return :(muladd( $(ks[i+1]*_f(x)), $((_a)(x, n, ks, i+1)) , 1))

# 2. Macro function as a for loop
macro _a(x,n,ks=k(n,α))
    ex = 1
    for i in n-1 : -1 : 0
        ex = :(muladd( $(ks[i+1]*_f(x)), $ex, 1))
    end
    return :($ex)
end

# 2. Second design, now with dispatching on types.
# 2.1 FIRST, we will define all the infrastructure for dispatching for Laguerre.
# 2.2 THEN, we will extend that design by adding dispatches for the Hermite case.
# 2.3 LAST, we will extend for all the other cases.

# 2.1 FIRST, the Laguerre case.
abstract type OP end
struct Laguerre{α,n} <: OP end
params(a :: T) where T = Tuple(T.parameters) # Credit to @yingboma
Laguerre(α, n :: Int) = Laguerre{α, n}()
Laguerre(n :: Int) = n >= 0 ? Laguerre{0,n}() : throw("The degree n must be non-negative.")
# Generalized binomial, useful too in the Jacobi case.
function Base.binomial(α :: T,k :: Int) where T<:AbstractFloat
    res = 1
    for i in 0:k-1
        res = res * (α - i)
    end
    return res/factorial(k)
end

d(p :: T, x) where T<:Laguerre = d(p)
function d(p :: T) where T<:Laguerre
    α,n = params(p)
    return binomial(α + n, n)
end

function b(p :: T, m) where T<:Laguerre
    α,n = params(p)
    return n-m+1
end

function c(p :: T, m) where T<:Laguerre
    α, n = params(p)
    return m*(α + m)
end

f(p :: T, x) where T<:Laguerre = x

# This now works for all Orthogonal Polynomials!
function k(p :: T, x) where T<:OP
    n = params(p)[end]
    return [-b(p, i)*inv(c(p, i)) for i in 1:n]
end

# Notice this works for all types Orthogonal Polynomials!
# Note: It doesn't return a specialized function yet - try to run it with @code_llvm
# We will make this specialized in the Third stage of design, after we add all the dispatches
# for all the polynomials.
function a(p :: T, x, m = params(p)[end]) where T<:OP
    res = oneunit(x)
    ks = k(p, x)

    for i in m-1: -1 : 0
        res = muladd(res, ks[i+1]*f(p, x),1)
    end

    return res * d(p, x)
end

# function a(p :: T, x, i=

# 2.2 SECOND - add a little abstract type machinery to handle dispatches for even/odd cases.
abstract type EvenOP <: OP end
abstract type OddOP <: OP end
struct Hermite{n} <: OP end
struct HermiteEven{n} <: EvenOP end
struct HermiteOdd{n} <: OddOP end

Hermite{n}() where n = iseven(n) ? HermiteEven{n ÷ 2}() : HermiteOdd{(n-1) ÷ 2}()
Hermite(n) =  n>= 0 ? Hermite{n}() : throw("The degree n must be non-negative.")

#Heven
function d(p :: T, x) where T<:HermiteEven
    n = params(p)[1]
    return (-1) ^ n * (factorial(2n)/factorial(n))
end

function b(p :: T, m) where T<:HermiteEven
    n = params(p)[1]
    return  2*(n - m + 1)
end

# This now covers half the cases!
c(p :: T, m) where T<:EvenOP = m*(2m - 1)

# This catches almost every case!
f(p :: T, x) where T<:OP = x^2

#Hodd
function d(p :: T, x) where T<:HermiteOdd
    n = params(p)[1]
    return (-1) ^ n * (factorial(2n +1)/factorial(n))*2x
end

function b(p :: T, m) where T<:HermiteOdd
    n = params(p)[1]
    return  2*(n - m + 1)
end

# This catches almost every case left over! We are only missing the Jacobi case.
c(p :: T, m) where T<:OddOP = m * (2m + 1)


# 2.3 THIRD, we now add dispatches to all the other cases.
struct Legendre{n} <: OP end
struct LegendreEven{n} <: EvenOP end
struct LegendreOdd{n} <: OddOP end
Legendre{n}() where n = iseven(n) ? LegendreEven{n ÷ 2}() : LegendreOdd{(n-1) ÷ 2}()
Legendre(n :: Int) = n>=0 ? Legendre{n}() : throw("The degree n must be non-negative.")

# PEven
# Note! We only have to define two functions per case here, because the other cases already got caught by
# some adequate planning and dispatching.
function d(p :: T, x) where T<:LegendreEven
    n = params(p)[1]
    return (-1/4)^n * binomial(2n,n)
end

function b(p :: T, m) where T<:LegendreEven
    n = params(p)[1]
    return (n - m + 1) * (2n + 2m -1)
end

# POdd
function d(p :: T, x) where T<:LegendreOdd
    n = params(p)[1]
    return (-1/4)^n * binomial(2n+1,n) * (n+1) * x
end

function b(p :: T, m) where T<:LegendreOdd
    n = params(p)[1]
    return (n - m + 1) * (2n + 2m + 1)
end


# TEven
struct ChebyshevFirstKind{n} <: OP end
struct ChebyshevFirstKindEven{n} <: EvenOP end
struct ChebyshevFirstKindOdd{n} <: OddOP end
ChebyshevFirstKind{n}() where n = iseven(n) ? ChebyshevFirstKindEven{n ÷ 2}() : ChebyshevFirstKindOdd{(n-1) ÷ 2}()
ChebyshevFirstKind(n :: Int) = n>=0 ? ChebyshevFirstKind{n}() : throw("The degree n must be non-negative")

function d(p :: T, x) where T<:ChebyshevFirstKindEven
    n = params(p)[1]
    return (-1)^n
end

function b(p :: T, m) where T<:ChebyshevFirstKindEven
    n = params(p)[1]
    return 2(n-m+1)*(n+m-1)
end


# TOdd
function d(p :: T, x) where T<: ChebyshevFirstKindOdd
    n = params(p)[1]
    return (-1^n)*(2n+1)*x
end

function b(p ::T, m) where T<: ChebyshevFirstKindOdd
    n = params(p)[1]
    return  2*(n-m+1)*(n+m)
end


# UEven
struct ChebyshevSecondKind{n} <: OP end
struct ChebyshevSecondKindEven{n} <: EvenOP end
struct ChebyshevSecondKindOdd{n} <: OddOP end
ChebyshevSecondKind{n}() where n = iseven(n) ? ChebyshevSecondKindEven{n ÷ 2}() : ChebyshevSecondKindOdd{(n-1) ÷ 2}()
ChebyshevSecondKind(n :: Int) = n>=0 ? ChebyshevSecondKind{n}() : throw("The degree n must be non-negative")

function d(p :: T, x) where T<:ChebyshevSecondKindEven
    n = params(p)[1]
    return (-1)^n
end

function b(p :: T, m) where T<:ChebyshevSecondKindEven
    n = params(p)[1]
    return 2*(n-m+1)*(n+m)
end

# UOdd
function d(p :: T, x) where T<:ChebyshevSecondKindOdd
    n = params(p)[1]
    return (-1)^n * 2(n+1)*x
end

function b(p :: T, m) where T<:ChebyshevSecondKindOdd
    n = params(p)[1]
    return 2*(n-m+1)*(n+m)
end


# CEven
struct Gegenbauer{α, n} <: OP end
struct GegenbauerEven{α, n} <: EvenOP end
struct GegenbauerOdd{α, n} <: OddOP end
function Gegenbauer(α, n :: Int)
     n >= 0 && α >= 0 ? nothing : throw("The degree n must be non-negative")
     iseven(n) ? GegenbauerEven{α, n ÷ 2}() : GegenbauerOdd{α, (n-1) ÷ 2}()
end
Gegenbauer(n :: Int) = Gegenbauer(0, n)
Gegenbauer

# We need to define a pochhammer symbol here!
function poch(a,::Val{N}) where N
            res = a; for i in Base.OneTo(N-1)
                res = res * (a+i)
            end
        return res
end
poch(a, ::Val{0}) = one(a)
poch(a, ::Val{1}) = a

function d(p :: T, x = 0) where T<:GegenbauerEven
    α, n  = params(p)
    return (-1)^n * (poch(α,Val(n))/factorial(n))
end

function b(p :: T, m) where T<:GegenbauerEven
    α, n = params(p)
    return 2*(n-m+1)*(α+n+m-1)
end

# COdd
function d(p :: T, x) where T<:GegenbauerOdd
    α, n = params(p)
    return (-1)^n * (poch(α,Val(n+1))/factorial(n))*2x
end

function b(p :: T, m) where T<:GegenbauerOdd
    α, n = params(p)
    return  2*(n-m+1)*(α + n + m)
end


# Jacobi
struct Jacobi{α, β, n} <: OP end
Jacobi(α, β, n :: Int) = α >= 0 && β >= 0 && n >= 0 ? Jacobi{α, β, n}() : throw("The parameters must be non-negative")

function d(p :: T, x) where T<:Jacobi
    α, β, n = params(p)
    return binomial(n + α, n)
end

function b(p :: T, m) where T<:Jacobi
    α, β, n = params(p)
    return (n - m + 1) * (α + β + n + m)
end

function c(p :: T, m) where T<:Jacobi
    α, β, n = params(p)
    return 2m * (α + m)
end

function f(p :: T, x) where T<:Jacobi
    return 1 - x
end

# Pass: Laguerre, ChebyFE, ChebySE,
# Check Jacobi-binomial, HermiteEven, HermiteOdd, LegendreE, LegendreO, ChebySO, ChebyFO, GegenbauerE, GegenbauerO
# Strat, look at the a(j,2,0) and k(j,.5) methods
export d,b,c,k,f,a
export Jacobi                                                                   # J(α,β,n)
export GegenbauerEven, GegenbauerOdd, Gegenbauer                                # C(α,n)
export ChebyshevSecondKindEven, ChebyshevSecondKindOdd, ChebyshevSecondKind     # U(n)
export ChebyshevFirstKindEven, ChebyshevFirstKindOdd, ChebyshevFirstKind        # T(n)
export Legendre                                                                 # L(n)
export HermiteOdd, HermiteEven, Hermite                                         # H(n)
export Laguerre

end # module
