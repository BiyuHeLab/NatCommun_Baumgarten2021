function o = nchoosek_g(n,k)

o = gamma(n+1) ./ ( gamma(k+1).*gamma(n-k+1) );