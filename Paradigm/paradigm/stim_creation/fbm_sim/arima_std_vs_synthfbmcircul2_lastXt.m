clear 

N      = 2^8+1;
betas  = linspace(1+1e-2,3-1e-2,10);
beta   = betas(1);
sigmas = [.1 1 3 5];
k      = 50;

nreps = 50;

for i_rep2 = 1:50
i_rep2
    for i_s = 1:length(sigmas)
   
    H_fBm = (beta-1)/2;
    sigma2 = sigmas(i_s)^2;

    for i_rep = 1:nreps
        [B x]  = synthfbmcircul2(N, H_fBm, sigma2);
        j = length(x);
        xe(i_rep) = x(j) - arima_pred(x(j-1:-1:j-k),H_fBm-.5);
        xe_t1(i_rep) = x(j-1) - arima_pred(x(j-1-1:-1:j-k-1),H_fBm-.5);
        
        xt1(i_rep) = x(j-1);
    end

    mxe(i_s) = mean(xe);
    sxe(i_s) = std(xe);

    r(i_rep2,i_s) = corr(xe',xt1');
    re(i_rep2,i_s) = corr(xe',xe_t1');
end

p=polyfit(sigmas,sxe,1);
m_farima(i_rep2) = p(1);

end

mean(m_farima)

fs = 15; 

figure; hold on;
plot(sigmas,sxe,'k-','LineWidth',2);
plot(sigmas,sigmas,'r-','LineWidth',2);

xlabel('generating sigma','FontSize',fs)
ylabel('recovered sigma (fARIMA, k = 50)','FontSize',fs)

axis([0 5 0 5]);
set(gca,'FontSize',fs)

