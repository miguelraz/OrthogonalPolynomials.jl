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
