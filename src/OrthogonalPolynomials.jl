__precompile__()

module OrthogonalPolynomials

# 1. First design - Laguerre only
d(n, α=0) = binomial(n+α, n)
b(n, m) = n-m+1
c(m, α) = m*(α+m)
k(n, α=0) = [-b(n,i)*inv(c(i, α)) * d(n, α) for i in 1:n]
f(x) = x
const α = 0

a(x, n, ks = k(n,α), i=0) = i == n ? :(1) : return :(muladd( $(ks[i+1]*f(x)), $((a)(x, n, ks, i+1)) , 1))

# 2. Macro function as a for loop
macro a(x,n,ks=k(n,α))
    ex = 1
    for i in n-1 : -1 : 0
        ex = :(muladd( $(ks[i+1]*f(x)), $ex, 1))
    end
    return :($ex)
end

export d,b,c,k,f,a,α

end # module
