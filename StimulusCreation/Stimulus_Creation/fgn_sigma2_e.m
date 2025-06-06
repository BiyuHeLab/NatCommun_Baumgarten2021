function sigma2_e = fgn_sigma2_e(H,sigma2,k)
% sigma2_e = fgn_sigma2_e(H,sigma2,k)
% 
% Calculate the theoretical sigma^2_{epsilon} from an fGn series with Hurst
% parameter H and overall variance sigma2, based on an estimate of the n'th
% sample of the fGn using the previous k samples

n = 0:k ;
r = sigma2/2*(abs(n-1).^(2*H) + abs(n+1).^(2*H) - 2*abs(n).^(2*H)) ;
Rbar = r(2:end)' ; 
temp = r(1:end-1) ; 
Rbarbar = toeplitz(temp') ; 
Abar = inv(Rbarbar)*Rbar  ;

sigma2_e = r(1) - sum(Abar'.*Rbar');