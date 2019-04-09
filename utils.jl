#Laguerre breaks at 16
for i in 12:20
    println("n = $i ", all(L_rec(i) .< 1/eps()))
end

# Laguerre Coefficient Calculator
# Credit to [Peter Luschny Apr. 11, 2015](https://oeis.org/A021009)
function l(n,m)
    if 0 == n == m
        1
    elseif m == -1
        0
    elseif n < m
        0
    elseif n >= 1
        (n+m+1) * l(n-1, m) - l(n-1, m-1)
    end
end

LTable(n) = [l(i,j) for i in 0:n,j in 0:n]

horner_L12(x) = @evalpoly x 1.0 -12.0 33.0 -36.666666666666664 20.625 -6.6 1.2833333333333334 -0.15714285714285714 0.012276785714285714 -0.0006062610229276896 1.8187830687830687e-5 -3.0062530062530064e-7

# Hermite Coefficient Calculator
# Credit to [Paul Barry, Aug 28, 2005.](https://oeis.org/search?q=coefficient+triangle+hermite&sort=&language=&go=Search)
function h(n,k)

    if iseven(n-k) && (n-k >= 0)
        ( (-1)^((n-k)/2) * (2^k) * factorial(n)) / (factorial(k) * factorial((n-k)/2))
    else
        0
    end
end

Htable(n) = [h(i,j) for i in 0:n, j in 0:n]

horner_H12(x) = Base.Math.@evalpoly x 665280 0 -7983360 0 13305600 0 -7096320 0 1520640 0 -135168 0 4096

# Legendre Coefficient Calculator
# Credit to Ralf Stephan, Apr. 07. 2016 https://oeis.org/A008556
PTable(n) =  [binomial(2*(i-j),i-j) * binomial(i-j,j) for i in 0:n, j in 0:n]

horner_P12(x) = Base.Math.@evalpoly x

# Tshebyshev T First Kind Coefficient Calculator
# Credit to Micahel Somos, Aug. 08, 2011 https://oeis.org/A049310
function T(n,k)
    if k < 0 || k > n || (n + k) % 2 == 1
        return 0
    end
    return (-1)^((n+k)/2 +k) * binomial(Int((n+k)/2),k)
end

TTable(n) = Int[u(i,j) for i in 0:n, j in 0:n]

# Tshebyshev U Second Kind Coefficient Calculator

function u(n,m)
    if n < m || isodd(n+m)
        0
    else
        ((-1)^((n+m)/2+m))*(2^m)*binomial(Int((n+m)/2), m)
    end
end

UTable(n) = Int[u(i,j) for i in 0:n, j in 0:n]

# Consider adding the Bell and Bernoulli polynomials


# stirling numbers of the first Kind

s1(m,n) = sum([(-1)^k*binomial(n-1+k,n-m+k)*binomial(2*n-m,n-m-k)*s2(k,n-m-k) for k in 0:n-m])
row_s1(n) = [s1(i,n) for i in 1:n]

# stirling numbers of the second kind

s2(n,k) = (1/factorial(k))*sum([(-1)^i * binomial(k,i) * (k - i)^n for i in 0:k])

row_s2(n) = [s2(n,i) for i in 1:n]

s(n,m) = sum([(-1)^(n-j) * binomial(k,j) * factorial(j*y - 1 + n) / factorial(j*y-1) * y ^ (-k) / factorial(k) for j in 0:k])

function stirlings1(n::Int, k::Int, signed::Bool=false)
    if signed == true
        return (-1)^(n - k) * stirlings1(n, k)
    end

    if n < 0
        throw(DomainError(n, "n must be nonnegative"))
    elseif n == k == 0
        return 1
    elseif n == 0 || k == 0
        return 0
    elseif n == k
        return 1
    elseif k == 1
        return factorial(n-1)
    elseif k == n - 1
        return binomial(n, 2)
    elseif k == n - 2
        return div((3 * n - 1) * binomial(n, 3), 4)
    elseif k == n - 3
        return binomial(n, 2) * binomial(n, 4)
    end

    return (n - 1) * stirlings1(n - 1, k) + stirlings1(n - 1, k - 1)
end

function stirlings2(n::Int, k::Int)
    if n < 0
        throw(DomainError(n, "n must be nonnegative"))
    elseif n == k == 0
        return 1
    elseif n == 0 || k == 0
        return 0
    elseif k == n - 1
        return binomial(n, 2)
    elseif k == 2
        return 2^(n-1) - 1
    end

    return k * stirlings2(n - 1, k) + stirlings2(n - 1, k - 1)
end

# Straight from Julia codebase
macro evalpoly(z, p...)
    a = :($(esc(p[end])))
    b = :($(esc(p[end-1])))
    as = []
    for i = length(p)-2:-1:1
        ai = Symbol("a", i)
        push!(as, :($ai = $a))
        a = :(muladd(r, $ai, $b))
        b = :($(esc(p[i])) - s * $ai) # see issue #15985 on fused mul-subtract
    end
    ai = :a0
    push!(as, :($ai = $a))
    C = Expr(:block,
             :(x = real(tt)),
             :(y = imag(tt)),
             :(r = x + x),
             :(s = muladd(x, x, y*y)),
             as...,
             :(muladd($ai, tt, $b)))
    R = Expr(:macrocall, Symbol("@horner"), (), :tt, map(esc, p)...)
    :(let tt = $(esc(z))
          isa(tt, Complex) ? $C : $R
      end)
end

#P
d(p, α,n) = binomial(n + α,n)
b(p, α, β, n, m) = (n-m+1)*(α+β+n+m)
c(p, α, m) = 2m*(α + m)
f(p, x) = 1 - x

#Ceven
d(p, α, n) = (-1)^n * (poch(α,n)/factorial(n))
b(p, α, n, m) = 2*(n-m+1)*(α+n+m-1)
c(p, m) = m*(2m-1)
f(p, x) = x^2

#Codd
d(p, α, n) = (-1)^n * (poch(α,n+1)/factorial(n))*2x
b(p, α, β, n, m) = 2*(n-m+1)*(α+β+n+m)
c(p, m) = m*(2m+1)
f(p, x) = x^2

#Teven
d(p, n) = (-1)^n
b(p, n, m) = 2(n-m+1)*(n+m-1)
c(p, m) = m*(2m-1)
f(p, x) = x^2

#Todd
d(p, n, x) = (-1^n)*(2n+1)*x
b(p, n, m) = 2*(n-m+1)*(n+m)
c(p, m) = m*(2m+1)
f(p, x) = x^2

#Ueven
d(p, n) = (-1)^n
b(p, n, m) = 2*(n-m+1)*(n+m)
c(p, m) = m*(2m+1)
f(p, x) = x^2

#Uodd
d(p, n, x) = (-1)^n * 2(n+1)*x
b(p, n, m) = 2*(n-m+1)*(n+m)
c(p, m) = m*(2m + 1)
f(p, x) = x^2

#Peven
d(p, n) = (-1/4)^n * binomial(2n,n)
b(p, n, m) = (n - m + 1) * (2n + 2m -1)
c(p, m) = m*(2m-1)
f(x) = x^2

#Podd
d(p, n, x) = (-1/4)^n * binomial(2n+1,n) * (n+1) * x
b(p, n, m) = (n - m + 1) * (2n + 2m + 1)
c(p, m) = m * (2m + 1)
f(p, x) = x^2

#L
d(p, α, n) = binomial(n + α, n)
b(p, n, m) = (n - m + 1)
c(p, m, α) = m * (α + m)
f(p, x) = x

#Heven
d(p, n) = (-1) ^ n * (factorial(2n)/factorial(n))
b(p, n, m) = 2*(n - m + 1)
c(p, m) = m*(2m - 1)
f(p, x) = x^2

#Hodd
d(p, n, x) = (-1) ^ n * (factorial(2n +1)/factorial(n))*2x
b(p, n, m) = 2*(n - m + 1)
c(p, m) = m * (2m + 1)
f(p, x) = x^2

#interpolate an expression as an expression into another expression.
ex_new = Meta.quot(ex)
quote
    still_expression = $(esc(ex_new))
end
