function [Xt Bkj] = arima_pred(Xts,d)
% [Xt Bkj] = arima_pred(Xts,d)
% 
% Estimate Xt from Xt-1, ..., Xt-k for an ARIMA(0,d,0) process, where
% input Xts is a vector of form [Xt-1, Xt-2, ..., Xt-k]
%
% Xt = sum(Bkj .* Xts)
%
% Based on Beran (1994) Statistics for Long Memory Processes, p 65

k = length(Xts);
for j=1:k
    Bkj(j) = -nchoosek(k,j) * (gamma(j-d) * gamma(k-d-j+1)) / (gamma(-d) * gamma(k-d+1));
end

Xt = sum(Bkj .* Xts);