function [Xt Bkj] = arima_pred2(Xts,d)
% [Xt Bkj] = arima_pred2(Xts,k,d)
% 
% Estimate Xt from Xt-1, ... Xt-k for an ARIMA(0,d,0) process, where
% input Xts is a vector of form [Xt-1, Xt-2, ..., Xt-k]
%
% Xt = sum(Bkj .* Xts)
%
% Based on Beran (1994) Statistics for Long Memory Processes, p 65

k = length(Xts);
ks = 1:k;

Fb = gamma(d + 1) ./ (gamma(ks + 1) .* gamma(d - ks + 1));
w  = -ones(size(Fb));
w  = w .^ ks;

Xt = -sum(w.*Fb.*Xts);