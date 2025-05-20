clear 

N     = 2^7+1;
betas = linspace(1+1e-2,3-1e-2,10);
k     = 20;

for i_b = 1:length(betas)
    
beta  = betas(i_b);
H_fBm = (beta-1)/2;
[B x]  = synthfbmcircul(N, H_fBm);

ind = 1;
for j = round(N/2) : N
    xe(ind) = x(j) - arima_pred(x(j-1:-1:j-k),H_fBm-.5);
    ind = ind+1; 
end

mxe(i_b) = mean(xe);
sxe(i_b) = std(xe);

end

H_fBms = (betas-1)/2;

figure; hold on;
plot(betas,sxe,'bo-');
% plot(betas,N.^-H_fBms,'ro-')
xlabel('generating beta for synthfbmcircul');
ylabel(['estimated sigma_e using k=' num2str(k)]);
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log'); 


% subplot(313); hold on;
figure; hold on;
plot(H_fBms,sxe,'bo-');
plot(H_fBms,N.^-H_fBms,'ro-')
xlabel('generating H_{fBm} for synthfbmcircul');
legend(['estimated sigma_e using k=' num2str(k)],'N^{-H}','Location','SouthWest');
title('N = 2^{7}+1')
set(gca,'XScale','log')
set(gca,'YScale','log')

