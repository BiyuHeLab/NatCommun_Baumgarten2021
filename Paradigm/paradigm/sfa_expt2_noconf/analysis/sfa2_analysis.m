clear 

% path = 'C:\Users\localbiyu\Desktop\sfa_expt2\data\original\';
path = 'C:\Users\maniscalcobs\Desktop\sfa_expt2\data\original\';

subs = {'BM' 'BH' 'KM' 'LS' 'KS' 'FC' 'JH' 'NM' 'CD' 'MP' 'LR' 'BS'};

sub = subs{8};
sub = 'BM';

switch sub
    case 'BM1', load([path 'Brian pilot 1\sfa_expt2_pilot brian.mat']);
    case 'BM', load([path 'Brian pilot 2\sfa_expt2 Brian pilot 2.mat']);
    case 'BH', load([path 'Biyu pilot 2\sfa_expt2 Biyu pilot 2.mat']);
    case 'KM', load([path 'KM\behavioral\sfa_expt2 2013-12-09 14-23-18.mat']);
    case 'LS', load([path 'LS\behavioral\sfa_expt2 2014-01-24 13-29-38.mat']);
    case 'KS', load([path 'KS\behavioral\sfa_expt2 2014-01-28 11-23-50.mat']);
    case 'FC', load([path 'FC\behavioral\sfa_expt2 2014-01-29 10-03-03.mat']);
    case 'JH', load([path 'JH\behavioral\sfa_expt2 2014-01-30 10-28-33.mat']);
    case 'NM', load([path 'NM\behavioral\sfa_expt2 2014-01-31 10-51-11.mat']);
    case 'CD', load([path 'CD\behavioral\sfa_expt2 2014-01-31 13-29-58.mat']);
    case 'MP', load([path 'MP\behavioral\sfa_expt2 2014-02-03 10-23-30.mat']);
    case 'LR', load([path 'LR\behavioral\sfa_expt2 2014-02-05 10-26-21.mat']);
    case 'BS', load([path 'BS\behavioral\sfa_expt2 2014-02-07 09-29-01.mat']);        
end


%% filter the data

manualFilter = ones(1,360);

switch sub
    case 'LS'
        % manual filter for LS
        % - omit trials in block 2 where r1_RT < .1 (button got stuck)
        % - omit trials in block 3 (button was stuck for most of block)
        manualFilter = manualFilter & data.r1_RT > .1;
        manualFilter = manualFilter & (data.trialNum < 61 | data.trialNum > 90);
end

% filter = data.resp_beta >= 0 & data.resp_prob >= 0 & data.conf_beta >= 0;
filter = data.resp_beta >= 0 & data.resp_prob >= 0 & data.conf_beta >= 0 & data.trialNum > 0;
filter = filter & data.r1_RT > .1 & data.r2_RT > .1 & data.r3_RT > .1;
filter = filter & manualFilter;

% filter = filter & stim.sigma_e_ID==0;
% filter = filter & abs(data.diff_beta) > .6;

dfn = fieldnames(data);
for j = 1:length(dfn)
    if eval(['length(data.' dfn{j} ') == 360'])
        eval(['dataf.' dfn{j} ' = data.' dfn{j} '(filter);']);
    end
end

sfn = fieldnames(stim);
for j = 1:length(sfn)
    if eval(['length(stim.' sfn{j} ') == 360'])
        eval(['stimf.' sfn{j} ' = stim.' sfn{j} '(filter);']);
    end
end


%% model selection

del_logf_pen  = abs( stimf.logf_final - stimf.logf_pen );
del_logf_pred = abs( stimf.logf_final - stimf.logf_pred );
del_SD_pen    = abs( stimf.logf_final - stimf.logf_pen ) ./ stimf.sigma_e;
del_SD_pred   = abs( stimf.logf_final - stimf.logf_pred ) ./ stimf.sigma_e;

dv = 5-dataf.resp_prob;
iv_list = {del_logf_pen,  del_logf_pred,  del_SD_pen,  del_SD_pred};
iv_labels2 = {'logf pen', 'logf pred', 'SD pen', 'SD pred'};
iv_labels = {'|log f_{t+1} - log f_t|', '|log f_{t+1} - E[log f_{t+1}]|', ...
             '|log f_{t+1} - log f_t| / \sigma_\epsilon', '|log f_{t+1} - E[log f_{t+1}]| / \sigma_\epsilon'};
% iv_labels2 = {'|log f_{t+1} - log f_t|', '|log f_{t+1} - E[log f_{t+1}|', ...
%              '|log f_{t+1} - log f_t| / s_e', '|log f_{t+1} - E[log f_{t+1}]| / s_e'};

         
h=figure;
fs = 15;
set(h,'DefaultAxesFontSize',fs);
for i_iv = 1:4
    iv = iv_list{i_iv};

    subplot(2,2,i_iv);
    plot(iv, dv, 'b.');
    ylim([.5 5.5])
    xlabel(iv_labels{i_iv});
    ylabel('prob rating')

end


for i_iv = 1:4
    iv = iv_list{i_iv};
    
    m = [dv' iv'];
    in.DATA = m;
    out = MATLAB_Ordered_Probit_Estimate(in);

    logL(i_iv) = out.Likelihood.LLV;

    % recover model probabilities
    x = iv*out.Beta;

    probit_x{i_iv} = x;
    probit_beta{i_iv} = out.Beta;
    probit_cutpoints{i_iv} = out.Cut_Points;

    for i_trial = 1:length(x)
        prob_model{i_iv}(i_trial) = 1 + sum( x(i_trial) > out.Cut_Points );
    end

% % %     [beta mu l] = fit_ordered_probit_MLE(iv', dv', 4);
% % %     logL(i_iv) = l;
% % %     x = iv*beta;
% % % 
% % %     probit_x{i_iv} = x;
% % %     probit_beta{i_iv} = beta;
% % %     probit_cutpoints{i_iv} = mu;
% % %     
% % %     for i_trial = 1:length(x)
% % %         prob_model{i_iv}(i_trial) = 1 + sum( x(i_trial) > mu );
% % %     end

   
% % %     [mu sd l] = fit_ordered_probit_MLE2(iv', dv', 4);
% % %     logL(i_iv) = l;
% % %     x = iv;
% % % 
% % %     probit_x{i_iv} = x;
% % %     probit_sd{i_iv} = sd;
% % %     probit_cutpoints{i_iv} = mu;
% % %     
% % %     for i_trial = 1:length(x)
% % %         prob_model{i_iv}(i_trial) = 1 + sum( x(i_trial) > mu );
% % %     end
    
    
    for k = 1:5
        p_b(k) = sum(dv==k)/length(dv);
        p_m{i_iv}(k) = sum(prob_model{i_iv}==k)/length(prob_model{i_iv});
    end
end

figure;
subplot(211);
b = [p_b' p_m{1}' p_m{2}' p_m{3}' p_m{4}'];
bar(b);

subplot(212);
bdiff = [p_b'-p_m{1}' p_b'-p_m{2}' p_b'-p_m{3}' p_b'-p_m{4}'];
bar(bdiff);

K = 5;
N = length(dv);

BIC = -2*logL + K.*log(N);
AIC = -2*logL + 2.*K.*N./(N-K-1);

BICw = exp( -.5*(BIC-min(BIC)) ) / sum( exp( -.5*(BIC-min(BIC)) ));
AICw = exp( -.5*(AIC-min(AIC)) ) / sum( exp( -.5*(AIC-min(AIC)) ))

figure;
fs = 15;

bar(AICw);
set(gca, 'XTickLabel', iv_labels2);
ylabel('Akaike weight', 'FontSize', fs)
set(gca, 'FontSize', fs);




%% beta discrimination

betas = unique(stimf.beta);
for i_b = 1:5
    for i_r = 1:5
        
        beta_stim = betas(i_b);
        beta_resp = betas(i_r);
        
        ff = stimf.beta == beta_stim;
        beta_matrix(i_b, i_r) = sum( dataf.resp_beta(ff) == beta_resp ) / length(dataf.resp_beta(ff));
 
        beta_cum(i_b, i_r) = sum( dataf.resp_beta(ff) <= beta_resp ) / length(dataf.resp_beta(ff));
        
        beta_nt(i_b, i_r) = sum(stimf.beta==beta_stim & dataf.resp_beta==beta_resp);
    end
end

for i_b = 2:5
    for i_r = 1:4
        beta_d(i_b-1, i_r) = norminv(beta_cum(i_b-1,i_r)) - norminv(beta_cum(i_b,i_r));
    end
end

[r p] = corr(stimf.beta', dataf.resp_beta', 'type', 'Spearman') 
[r1 p1] = corr(stimf.sigma_e', dataf.resp_beta', 'type', 'Spearman') 

% [r2 p2] = partialcorr(dataf.resp_beta', stimf.beta', stimf.sigma_e','type','Spearman')
% [r3 p3] = partialcorr(dataf.resp_beta', stimf.sigma_e', stimf.beta','type','Spearman')
% 
% ff = stimf.sigma_e_ID==1 & stimf.beta < 1.5;
% [r4 p4] = corr(stimf.beta(ff)', dataf.resp_beta(ff)', 'type', 'Spearman') 
% 
% ff = stimf.sigma_e_ID==0 & stimf.beta < 1.5;
% [r5 p5] = corr(stimf.beta(ff)', dataf.resp_beta(ff)', 'type', 'Spearman') 
% 
% ff = (stimf.sigma_e_ID==0 & stimf.beta==0) | (stimf.sigma_e_ID==1 & stimf.beta==2);
% [r5 p5] = corr(stimf.beta(ff)', dataf.resp_beta(ff)', 'type', 'Spearman') 



figure;
fs = 15;

imagesc(beta_matrix)
axis square
colormap hot
colorbar

xlabel('response beta', 'FontSize', fs)
ylabel('stimulus beta', 'FontSize', fs)
title(['p( response beta | stim beta ), rho = ' num2str(r,2)], 'FontSize', fs)

set(gca,'XTick',1:5)
set(gca,'YTick',1:5)
set(gca,'XTickLabel',{'0' '0.5' '1' '1.5' '2'})
set(gca,'YTickLabel',{'0' '0.5' '1' '1.5' '2'})
set(gca,'FontSize',fs)

caxis([0 .71])




figure;
fs = 15;

% imagesc(beta_matrix)
imagesc(beta_cum)
axis square
colormap hot
colorbar

xlabel('response beta', 'FontSize', fs)
ylabel('stimulus beta', 'FontSize', fs)
title(['p( response beta | stim beta ), rho = ' num2str(r,2)], 'FontSize', fs)

set(gca,'XTick',1:5)
set(gca,'YTick',1:5)
set(gca,'XTickLabel',{'0' '0.5' '1' '1.5' '2'})
set(gca,'YTickLabel',{'0' '0.5' '1' '1.5' '2'})
set(gca,'FontSize',fs)

% caxis([0 .71])
caxis([0 1])




figd; hold on;

beta_d(abs(beta_d)==Inf) = NaN;
dd = nanmean(beta_d,2);
dcum = cumsum(dd');

plot([0 .5:.5:2], [0 dcum], 'bo');
pp = polyfit(.5:.5:2, dcum, 1);
xx = [.5 2];
plot(xx, xx*pp(1)+pp(2), 'b-'); 
xlabel('beta')
ylabel('cumulative d''');
xlim([-.2 2.2])