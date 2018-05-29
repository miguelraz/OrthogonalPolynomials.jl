

#Laguerre breaks at 16
for i in 12:20
    println("n = $i ", all(L_rec(i) .< 1/eps()))
end



# Laguere Coefficient Calculator
# Credit to [Peter Luschny Apr. 11, 2015](https://oeis.org/search?q=laguerre+coefficient+triangle&sort=&language=&go=Search)
function r(n,m)
    if 0 == n == m
        1
    elseif m == -1
        0
    elseif n < m
        0
    elseif n >= 1
        (n+m+1) * r(n-1, m) - r(n-1, m-1)
    end
end

L_rec(n) = [r(j,i) for i in 0:n,j in 0:n]

# Credit to [Paul Barry, Aug 28, 2005.](https://oeis.org/search?q=coefficient+triangle+hermite&sort=&language=&go=Search)
# Hermite Coefficient Calculator
function T(n,k)

    if iseven(n-k) && (n-k >= 0)
        ( (-1)^((n-k)/2) * (2^k) * factorial(n)) / (factorial(k) * factorial((n-k)/2))
    else
        0
    end
end
