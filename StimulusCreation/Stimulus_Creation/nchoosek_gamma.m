function out = nchoosek_gamma(n,k)
% out = nchoosek_g(n,k)
% 
% Calculate nhoosek using the gamma function
%
% out = gamma(n+1) ./ ( gamma(k+1).*gamma(n-k+1) );

out = gamma(n+1) ./ ( gamma(k+1).*gamma(n-k+1) );