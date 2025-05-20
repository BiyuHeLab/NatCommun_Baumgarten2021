H = 1e-2;
d = H - .5;
N = 2^10+1;

[B x] = synthfbmcircul(N,H);

Xt  = x(end);

k_list = 1:55;

for i_k = 1:length(k_list)
    
Xts = x(end-1 : -1 : end - k_list(i_k));

Xt_est(i_k) = arima_pred(Xts,d);
% Xt_est2(i_k) = arima_pred2(Xts,d);
% Xt_est3(i_k) = arima_pred3(Xts,d);

end

figure; hold on;
plot(k_list, Xt-Xt_est,'b-');
% plot(k_list, Xt-Xt_est2,'r-');
% plot(k_list, Xt-Xt_est3,'g-');