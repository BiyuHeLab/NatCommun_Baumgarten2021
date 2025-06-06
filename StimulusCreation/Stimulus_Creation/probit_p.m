function prob = probit_p(x, beta, mu)
% prob = probit_p(x, beta, mu)
%
% Given an independent variable x, weighting factor beta, and thresholds mu, 
% calculate the probability of ordinal output variable 
%
% y = j, for j = 0 to m.
%
% where y* = x*beta + e, and y = j iff y* is located between the
% appropriate thresholds in mu.
%
% Inputs
% ------
% x    - 1 x V vector of independent variables
% beta - V x 1 vector of weights
% mu   - 1 x m vector of thresholds
%
% Output
% ------
% prob - 1 x m+1 vector of probabilities. 
%        For each entry i, prob(i) = p( y == i-1 | x, beta, mu )


mu = [-Inf mu Inf];
m  = length(mu) - 2;

prob = normcdf( mu(2:end) - x*beta ) - normcdf( mu(1:end-1) - x*beta );

end