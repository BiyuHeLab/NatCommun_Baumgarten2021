clear

nrep = 1e4;

H = 1-1e-2;
d = H - .5;
N = 2^10;

Xt = zeros(1,nrep);
% Xts = zeros(k,nrep);
% Bkj = zeros(size(Xts));

k_list = [1 3 5 10];

for i_k=1:length(k_list)

i_k
    
for i=1:nrep

[B x] = synthfbmcircul(N,H);

% Xt(i) = B(end);
% Xts   = B(end-1 : -1 : end-k);

Xt(i) = x(end);
Xts   = x(end-1 : -1 : end - k_list(i_k));

Xt_est(i) = arima_pred(Xts,d);
Xt_est2(i) = arima_pred2(Xts,d);

end

Xt_err = Xt - Xt_est;
Xte_mean(i_k) = mean(Xt_err);
Xte_std(i_k)  = std(Xt_err);

Xt_err2 = Xt - Xt_est2;
Xte_mean2(i_k) = mean(Xt_err2);
Xte_std2(i_k)  = std(Xt_err2);

end

figure;
subplot(211); hold on;
plot(k_list,Xte_mean,'bo-');
plot(k_list,Xte_mean2,'ro-');
xlabel('k');
ylabel('Xt error mean');

subplot(212); hold on;
plot(k_list,Xte_std,'bo-');
plot(k_list,Xte_std2,'ro-');
xlabel('k');
ylabel('Xt error std');



% hist(Xt_err);
% title(['mean = ' num2str(mean(Xt_err)) ', std = ' num2str(std(Xt_err)) ', k = ' num2str(k)]);
% 
% subplot(211);
% plot(B);
% subplot(212);
% plot(x);