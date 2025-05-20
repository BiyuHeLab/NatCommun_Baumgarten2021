clear

betas  = [0 .25 .5 .75 1.01 1.25 1.5 1.75 2 2.25 2.5 2.75 2.9];
betas = [0 : .1 : 2.9];
betas(betas==1) = 1.01; 

sigmas = [.5 1 5 10];

nSeries = 1000;
k = 50;

for i_b = 1:length(betas)
    t0 = GetSecs;
    percent = (i_b-1)*100/length(betas)
    for i_s = 1:length(sigmas)
        
        beta   = betas(i_b);
        sigma2 = sigmas(i_s)^2;
        [s sp pe fGns] = converging_prediction(beta, sigma2, k, 0, nSeries);
        
        if beta > 1, H = (beta-1)/2; else H = (beta+1)/2; end
        sigma2_e = fgn_sigma2_e(H,sigma2,k);
        
        sigma_pred(i_b,i_s) = std(sp);
        sigma_epsilon(i_b,i_s) = sqrt(sigma2_e);
        sigma_epsilon_sim(i_b,i_s) = std(pe);
        
        for i_f = 1:nSeries
            fGn = fGns(i_f,:);
            std_fGn(i_f) = std(fGn);
        end
        sigma_fGn(i_b,i_s) = mean(std_fGn);
        % sigma_fGn(i_b,i_s) / sqrt(sigma2)
        
    end
    tf(i_b) = (GetSecs-t0)/60
end


% prediction sigma vs fgn sigma
% betas  = [.25 .5 .75];
    
figure; hold on;
Rs = linspace(0, 1, length(betas));
Bs = 1-Rs;
for i_b = 1:length(betas)
    beta = betas(i_b)
    sigma_pred(i_b,:) ./ sigmas
    plot(sigmas.^2, sigma_pred(i_b,:).^2,'b^-','Color',[Rs(i_b) 0 0],'LineWidth',2);
%     plot(sigmas, sigma_epsilon(i_b,:),'bv-','Color',[Rs(i_b) 0 0],'LineWidth',2);

end
xlabel('fGn sigma');
% ylabel('std of fGn prediction for next sample')



% variance parcellation
for i_b = 1:length(betas)
    p = polyfit(sigmas.^2, sigma_pred(i_b,:).^2,1);
    m_s_pred(i_b) = p(1);
    
    p = polyfit(sigmas.^2, sigma_epsilon(i_b,:).^2,1);
    m_s_epsilon(i_b) = p(1);
    
    p = polyfit(sigmas.^2, sigma_epsilon_sim(i_b,:).^2,1);
    m_s_epsilon_sim(i_b) = p(1);
    
    p = polyfit(sigmas.^2, sigma_fGn(i_b,:).^2,1);
    m_s_fGn(i_b) = p(1);

end

figure; hold on;
plot(betas,m_s_pred,'bo-');
plot(betas,m_s_epsilon,'ro-');
plot(betas,m_s_pred+m_s_epsilon,'go-');

plot(betas,ones(size(betas)),'k-');
plot(betas,m_s_fGn,'ko-');

% plot(betas,1-m_s_epsilon,'bd--');
% plot([betas],[m_s_epsilon_sim],'r*--');
% plot([betas],[m_s_pred+m_s_epsilon_sim],'g*--');


xlabel('beta');
% ylabel('sigma of predicted (n+1)''th sample / sigma of fGn')