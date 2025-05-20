clear 

N      = 2^8+1;
betas  = linspace(1+1e-2, 3-1e-2, 10);
sigmas = [.1 1 3 5];
k      = 50;
k2 = 50;

nreps = 50;

for i_rep = 1:nreps
    i_rep
    for i_b = 1:length(betas)
        for i_s = 1:length(sigmas)

            beta   = betas(i_b);
            H_fBm  = (beta-1)/2;
            sigma2 = sigmas(i_s)^2;
            [B x]  = synthfbmcircul2(N, H_fBm, sigma2);
            
            ind = 1;
            for j = round(N/2) : N
                % xe = diff b/t observed x and farima-predicted x at time j
                % we collect many xe's and calculate the std of these to estimate
                % sigma_e as defined in Beran p 59
                xe(ind) = x(j) - arima_pred(x(j-1:-1:j-k), H_fBm-.5);
                
                [xhat A sigma2_e] = fgn_pred(x(j-1:-1:j-k2), H_fBm, sigma2);
                xe2(ind) = x(j) - xhat;
                
                ind = ind+1; 
            end

            mxe(i_b,i_s)   = mean(xe);
            sxe(i_b,i_s)   = std(xe);
            
            mxe2(i_b,i_s)   = mean(xe2);
            sxe2(i_b,i_s)   = std(xe2);
            
            std_x(i_b,i_s) = std(x);
            sig_e(i_b,i_s) = sqrt(sigma2_e);
        end

        p=polyfit(sigmas,sxe(i_b,:),1);
        m_farima(i_rep,i_b) = p(1);
        
        p=polyfit(sigmas,sxe2(i_b,:),1);
        m_fgn(i_rep,i_b) = p(1);        

        p=polyfit(sigmas,std_x(i_b,:),1);
        m_stdx(i_rep,i_b) = p(1);
        
        p=polyfit(sigmas,sig_e(i_b,:),1);
        m_sige(i_rep,i_b) = p(1);
        

    end
end



fs = 15; 

% plot recovered sigmas vs generating sigma for farima and fGn methods
% only plots data from the final simulation repetition
figure; 
subplot(1,3,1); hold on;
color = linspace(0,.8,10);
for i_b=1:length(betas)
    c = color(i_b);
    plot(sigmas,sxe(i_b,:),'Color',[c c c],'LineWidth',2);
end
plot(sigmas,sigmas,'r-','LineWidth',2);
xlabel('generating sigma','FontSize',fs)
ylabel(['recovered sigma (fARIMA, k = ' num2str(k) ')'],'FontSize',fs)
title('lighter lines = higher beta (ranging 1 - 3)','FontSize',fs)
axis([0 5 0 5]);
set(gca,'FontSize',fs)

subplot(1,3,2); hold on;
color = linspace(0,.8,10);
for i_b=1:length(betas)
    c = color(i_b);
    plot(sigmas,std_x(i_b,:),'Color',[c c c],'LineWidth',2);
end
plot(sigmas,sigmas,'r-','LineWidth',2);
xlabel('generating sigma','FontSize',fs)
ylabel('recovered sigma (std(x))','FontSize',fs);
axis([0 5 0 5]);
set(gca,'FontSize',fs)

subplot(1,3,3); hold on;
color = linspace(0,.8,10);
for i_b=1:length(betas)
    c = color(i_b);
    plot(sigmas,sxe2(i_b,:),'Color',[c c c],'LineWidth',2);
end
plot(sigmas,sigmas,'r-','LineWidth',2);
xlabel('generating sigma','FontSize',fs)
ylabel(['recovered sigma (fgn pred, k = ' num2str(k2) ')'],'FontSize',fs);
axis([0 5 0 5]);
set(gca,'FontSize',fs)


% plot average slopes (across simulation repetitions) for recovered vs 
% generating sigma as a function of beta
figure; hold on;
errorbar(betas,mean(m_stdx),std(m_stdx)/sqrt(nreps),'r-','LineWidth',2);
plot(betas,m_sige(1,:),'ko-','LineWidth',2);
errorbar(betas,mean(m_fgn),std(m_fgn)/sqrt(nreps),'g-','LineWidth',2);
errorbar(betas,mean(m_farima),std(m_farima)/sqrt(nreps),'b-','LineWidth',2);

xlabel('beta','FontSize',fs)
ylabel('mean regression slope for IV = sigma used to generate the fGn','FontSize',fs)
set(gca,'FontSize',fs)
legend('std(x)','fGn sigma_e (theoretical)',['fGn sigma_e (k=' num2str(k2) ')'],['fARIMA sigma_e(k=' num2str(k) ')'],'FontSize',fs,'Location','SouthWest')
title(['n=' num2str(nreps) ' repetitions'],'FontSize',fs);

