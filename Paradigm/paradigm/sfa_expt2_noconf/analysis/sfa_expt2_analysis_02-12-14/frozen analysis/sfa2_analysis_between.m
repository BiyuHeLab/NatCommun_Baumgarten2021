%% set up
clear 

model_selection = 1;
do_recov = 0;

% path = 'C:\Users\localbiyu\Desktop\sfa_expt2\data\original\';
path = 'C:\Users\maniscalcobs\Desktop\sfa_expt2\data\original\';

subs = {'BM' 'BH' 'KM' 'LS' 'KS' 'FC' 'JH' 'NM' 'CD' 'MP' 'LR' 'BS'};

% % remove pilot subjects BM and BH
% subs = subs(3:end);

%%% remove poorly performing subjects
% pilot subjects
bad_pilot = {'BM' ... % pilot
             'BH' ... % pilot
            }; 
   
bad_perf = {'JH' ... % maximum cumulative d' = 0.6
            'FC' ... % maximum cumulative d' = 0.3!
           };

bad_trial = {'LS' ... % lost 104 trials
             'NM' ... % lost 147 trials
            };


subs = setdiff(subs, bad_pilot);
subs = setdiff(subs, bad_perf);
% subs = setdiff(subs, bad_trial);

ss = 'BM1';
subs = {ss ss};

nSubs = length(subs);


%% loop over subjects

for i_sub = 1:length(subs)
    
sub = subs{i_sub};

switch sub
    case 'BM1'
        load([path 'Brian pilot 1\sfa_expt2_pilot brian.mat']);
        data.conf_beta = data.resp_prob;
        stim.betaID = ones(1,360);
        stim.betaID(stim.beta==.5) = 2;
        stim.betaID(stim.beta==1.01) = 3;
        stim.betaID(stim.beta==1.5) = 4;
        stim.betaID(stim.beta==2) = 5;
        data.r1_RT = data.r1_responseRT;
        data.r2_RT = data.r2_responseRT;
        data.r3_RT = data.r2_RT;
        
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

% figure; 
% hold on;
% plot(data.resp_prob, 'b-')
% plot(data.conf_beta, 'r-')
% title(sub)

%% adjust and define some data variables for later processing

% adjust diff_beta so that it only takes values in 0: .5 : 2
% (i.e. correct for subtractions by beta = 1.01)
data.diff_beta = round( 10*data.diff_beta ) / 10;

% adjust diff_beta so that it only takes values in 0: .5 : 2
% (i.e. correct for subtractions by beta = 1.01)
data.diff_beta3 = data.diff_beta;
data.diff_beta3(data.diff_beta3 < -1) = -1;
data.diff_beta3(data.diff_beta3 > +1) = +1;

% make an index for diff_beta3_ID
data.diff_beta3_ID( abs(data.diff_beta3) == 1 )  = 1;
data.diff_beta3_ID( abs(data.diff_beta3) ==.5 )  = 2;
data.diff_beta3_ID( abs(data.diff_beta3) == 0 )  = 3;

% make an index for resp_beta
data.resp_beta_ID(data.resp_beta==   0) = 1;
data.resp_beta_ID(data.resp_beta==  .5) = 2;
data.resp_beta_ID(data.resp_beta==1.01) = 3;
data.resp_beta_ID(data.resp_beta== 1.5) = 4;
data.resp_beta_ID(data.resp_beta==   2) = 5;


%% add in data for "recovered" stim parameters related to final tone estimation
load('series_tf_recov.mat');

% names of recov vars stored in series_tf_recov.mat
recov_vars = {'sigma_recov', 'H_recov', 'sigma_e_recov', 'tf_pred_recov', ...
    'sigma_recov_sc', 'H_recov_sc', 'sigma_e_recov_sc', 'tf_pred_recov_sc'};


% names to be used in stim.recov
recov_vars_newnames = {'sigma_recov', 'H_recov', 'sigma_e_recov', 'logf_pred_recov'};


for i_r = 1:length(recov_vars_newnames)
    eval(['stim.recov.' recov_vars_newnames{i_r} ' = NaN(1,360);']);
end

% for each trial, find the values for sigma, H, sigma_e, and the prediction
% for the final tone, based on what the "ideal observer" could recover 
% from the presented sequence...
for i_t = 1:360
    
    % get beta index for this trial
    i_b = stim.betaID(i_t);

    
    % get sequence index for this trial
    switch num2str(stim.predID(i_t))
        case '-1', i_s = 1; % min predicted frequency
        case '1',  i_s = 2; % max predicted frequency
        case '0',  i_s = 3; % median predicted frequency
    end
    
    
    % get sigma_e index for this trial,
    % which determines if we will use the first 4 "recov" variables (high sigma_e)
    % or the last 4 (low sigma_e)
    if stim.sigma_e_ID == 1
        i_r_offset = 0;
    else
        i_r_offset = 4;
    end
    
    for i_r  = 1:4
        eval(['stim.recov.' recov_vars_newnames{i_r} '(i_t) = ' recov_vars{i_r + i_r_offset} '(i_b, i_s);']);
    end
end
    

%% filter the data

filter_manual = ones(1,360);

switch sub
    case 'LS'
        % manual filter for LS
        % - omit trials in block 2 where r1_RT < .1 (button got stuck)
        % - omit trials in block 3 (button was stuck for most of block)
%         filter_manual = filter_manual & data.r1_RT > .1;
        filter_manual = filter_manual & (data.trialNum < 61 | data.trialNum > 90);

    case 'LR'
%         % LR appears to have swapped the final tone probability rating...
%         % tho maybe not... swapping her resp_prob doesn't improve her
%         % model comparison results at all...
%         data.resp_prob = 6 - data.resp_prob;
        
        
end

filter_resp  = data.resp_beta >= 0 & data.resp_prob >= 0 & data.conf_beta >= 0;
filter_trial = data.trialNum > 0;
filter_rt    = data.r1_RT > .1 & data.r2_RT > .1 & data.r3_RT > .1;
filter       = filter_resp & filter_trial & filter_rt & filter_manual;

nTrialsLost_manual(i_sub) = sum(filter_manual==0);
nTrialsLost_resp(i_sub)   = sum(filter_resp==0);
nTrialsLost_trial(i_sub)  = sum(filter_trial==0);
nTrialsLost_rt(i_sub)     = sum(filter_rt==0);
nTrialsLost(i_sub)        = sum(filter==0);

% filter the data struct
dfn = fieldnames(data);
for j = 1:length(dfn)
    if eval(['length(data.' dfn{j} ') == 360'])
        eval(['dataf.' dfn{j} ' = data.' dfn{j} '(filter);']);
    end
end


% filter the stim struct
sfn = fieldnames(stim);
for j = 1:length(sfn)
    if eval(['length(stim.' sfn{j} ') == 360'])
        eval(['stimf.' sfn{j} ' = stim.' sfn{j} '(filter);']);
    end
end


% filter the stim.recov struct
rfn = fieldnames(stim.recov);
for j = 1:length(rfn)
    if eval(['length(stim.recov.' rfn{j} ') == 360'])
        eval(['stimf.recov.' rfn{j} ' = stim.recov.' rfn{j} '(filter);']);
    end
end


%% beta discrimination

betas = unique(stimf.beta);
for i_b = 1:5
    for i_r = 1:5
        
        beta_stim = betas(i_b);
        beta_resp = betas(i_r);
        
        ff = stimf.beta == beta_stim;
        beta_matrix(i_b, i_r, i_sub) = sum( dataf.resp_beta(ff) == beta_resp ) / length(dataf.resp_beta(ff));
 
        beta_cum(i_b, i_r, i_sub) = sum( dataf.resp_beta(ff) <= beta_resp ) / length(dataf.resp_beta(ff));
        
        beta_nt(i_b, i_r, i_sub) = sum(stimf.beta==beta_stim & dataf.resp_beta==beta_resp);
    end
end

for i_b = 2:5
    for i_r = 1:4
        beta_d(i_b-1, i_r, i_sub) = norminv(beta_cum(i_b-1,i_r,i_sub)) - norminv(beta_cum(i_b,i_r,i_sub));
    end
end

[r_beta(i_sub) p_beta(i_sub)] = corr(stimf.beta', dataf.resp_beta', 'type', 'Spearman');

f = stimf.beta < 1.5;
[r_beta_low(i_sub) p_beta_low(i_sub)] = corr(stimf.beta(f)', dataf.resp_beta(f)', 'type', 'Spearman');



%% confidence rating

for i_b = 1:5
    % confidence for each stimulus beta
    f = stimf.betaID == i_b;
    conf(i_sub, i_b) = mean(dataf.conf_beta(f));
    
    % confidence vs beta for correct responses
    f = stimf.betaID == i_b & dataf.correct_beta == 1;
    conf_corr1(i_sub, i_b) = mean(dataf.conf_beta(f));

    % confidence vs beta for incorrect responses
    f = stimf.betaID == i_b & dataf.correct_beta == 0;
    conf_corr0(i_sub, i_b) = mean(dataf.conf_beta(f));
    
    
    % confidence vs beta for different levels of diff_beta     
    for i_diff = 1:3          
        f = stimf.betaID == i_b & dataf.diff_beta3_ID == i_diff;        
        conf_diff(i_diff, i_b, i_sub) = mean(dataf.conf_beta(f));
    end
    
    
    % conf for each stim/resp beta
    for i_r = 1:5
        f = stimf.betaID == i_b & dataf.resp_beta_ID == i_r;
        conf_br(i_r,i_b,i_sub) = mean(dataf.conf_beta(f));
    end
    
end



%% conf rating matrix

for i_conf = 1:5
    for i_diff = 1:3
              
        ff = dataf.diff_beta3_ID == i_diff;
        conf_matrix(i_diff, i_conf, i_sub) = sum( dataf.conf_beta(ff) == i_conf ) / length(dataf.conf_beta(ff));
 
        conf_nt(i_diff, i_conf, i_sub) = sum(dataf.conf_beta==i_conf & dataf.diff_beta3_ID==i_diff);
    end
end

[r_conf(i_sub)  p_conf(i_sub) ] = corr(dataf.conf_beta', abs(dataf.diff_beta'),  'type', 'Spearman');
[r_conf3(i_sub) p_conf3(i_sub)] = corr(dataf.conf_beta', abs(dataf.diff_beta3)', 'type', 'Spearman');



% repeat for each level of beta...
for i_b = 1:5
    for i_conf = 1:5
        for i_diff = 1:3

            ff = dataf.diff_beta3_ID == i_diff & stimf.betaID == i_b;
            conf_matrix_b{i_b}(i_diff, i_conf, i_sub) = sum( dataf.conf_beta(ff) == i_conf ) / length(dataf.conf_beta(ff));

            conf_nt_b{i_b}(i_diff, i_conf, i_sub) = sum(dataf.conf_beta==i_conf & dataf.diff_beta3_ID==i_diff);
        end
    end
    
    f = stimf.betaID == i_b;
    [r_conf_b{i_b}(i_sub)  p_conf_b{i_b}(i_sub) ] = corr(dataf.conf_beta(f)', abs(dataf.diff_beta(f))',  'type', 'Spearman');
    [r_conf3_b{i_b}(i_sub) p_conf3_b{i_b}(i_sub)] = corr(dataf.conf_beta(f)', abs(dataf.diff_beta3(f))', 'type', 'Spearman');

end




% repeat for each level of reported beta...
for i_r = 1:5
    for i_conf = 1:5
        for i_diff = 1:3

            ff = dataf.diff_beta3_ID == i_diff & dataf.resp_beta_ID == i_r;
            conf_matrix_r{i_r}(i_diff, i_conf, i_sub) = sum( dataf.conf_beta(ff) == i_conf ) / length(dataf.conf_beta(ff));

            conf_nt_r{i_r}(i_diff, i_conf, i_sub) = sum(dataf.conf_beta==i_conf & dataf.diff_beta3_ID==i_diff);
        end
    end
    
    f = dataf.resp_beta_ID == i_r;
    [r_conf_r{i_r}(i_sub)  p_conf_r{i_r}(i_sub) ] = corr(dataf.conf_beta(f)', abs(dataf.diff_beta(f))',  'type', 'Spearman');
    [r_conf3_r{i_r}(i_sub) p_conf3_r{i_r}(i_sub)] = corr(dataf.conf_beta(f)', abs(dataf.diff_beta3(f))', 'type', 'Spearman');

end

%% final tone prediction

% correlate final tone prob rating with the 4 permutations of predicted / penultimate, logf / sd

if ~do_recov
    % undo the recov...
    stimf.recov.logf_pred_recov = stimf.logf_pred;
    stimf.recov.sigma_e_recov   = stimf.sigma_e;
end

abs_pen  = abs( stimf.logf_final - stimf.logf_pen );
abs_pred = abs( stimf.logf_final - stimf.recov.logf_pred_recov );
sd_pen   = abs_pen  ./ stimf.recov.sigma_e_recov;
sd_pred  = abs_pred ./ stimf.recov.sigma_e_recov;

f  = dataf.diff_beta > -Inf; % arbitrary filter condition which selects all trials
% f  = abs(dataf.diff_beta) > .5;
% f = dataf.conf_beta <= median(dataf.conf_beta);

corr_type = 'Pearson';
[r_abs_pen(i_sub)  p_abs_pen(i_sub) ] = corr(dataf.resp_prob(f)', abs_pen(f)',  'type', corr_type);
[r_abs_pred(i_sub) p_abs_pred(i_sub)] = corr(dataf.resp_prob(f)', abs_pred(f)', 'type', corr_type);
[r_sd_pen(i_sub)   p_sd_pen(i_sub)  ] = corr(dataf.resp_prob(f)', sd_pen(f)',   'type', corr_type);
[r_sd_pred(i_sub)  p_sd_pred(i_sub) ] = corr(dataf.resp_prob(f)', sd_pred(f)',  'type', corr_type);


% do the same, for each level of beta...
for i_b = 1:5
    f_beta = stimf.betaID == i_b;
    
    abs_pen  = abs( stimf.logf_final - stimf.logf_pen );
    abs_pred = abs( stimf.logf_final - stimf.recov.logf_pred_recov );
    sd_pen   = abs_pen  ./ stimf.recov.sigma_e_recov;
    sd_pred  = abs_pred ./ stimf.recov.sigma_e_recov;

    abs_pen_f   = abs_pen(f_beta);
    abs_pred_f  = abs_pred(f_beta);
    sd_pen_f    = sd_pen(f_beta);
    sd_pred_f   = sd_pred(f_beta);
    resp_prob_f = dataf.resp_prob(f_beta);
    
    corr_type = 'Pearson';
    [r_abs_pen_b(i_sub, i_b)  p_abs_pen_b(i_sub, i_b) ] = corr(resp_prob_f', abs_pen_f',  'type', corr_type);
    [r_abs_pred_b(i_sub, i_b) p_abs_pred_b(i_sub, i_b)] = corr(resp_prob_f', abs_pred_f', 'type', corr_type);
    [r_sd_pen_b(i_sub, i_b)   p_sd_pen_b(i_sub, i_b)  ] = corr(resp_prob_f', sd_pen_f',   'type', corr_type);
    [r_sd_pred_b(i_sub, i_b)  p_sd_pred_b(i_sub, i_b) ] = corr(resp_prob_f', sd_pred_f',  'type', corr_type);    
end


%% model selection

iv_labels = {'|log f_{t+1} - log f_t|', '|log f_{t+1} - E[log f_{t+1}]|', ...
             '|log f_{t+1} - log f_t| / \sigma_\epsilon', '|log f_{t+1} - E[log f_{t+1}]| / \sigma_\epsilon'};


if model_selection

if ~do_recov
    % undo the recov...
    stimf.recov.logf_pred_recov = stimf.logf_pred;
    stimf.recov.sigma_e_recov   = stimf.sigma_e;
end

abs_pen  = abs( stimf.logf_final - stimf.logf_pen );
abs_pred = abs( stimf.logf_final - stimf.recov.logf_pred_recov );
sd_pen   = abs_pen  ./ stimf.recov.sigma_e_recov;
sd_pred  = abs_pred ./ stimf.recov.sigma_e_recov;

del_logf_pen  = abs_pen;
del_logf_pred = abs_pred;
del_SD_pen    = sd_pen;
del_SD_pred   = sd_pred;

% resp_prob ranges from 1-5 in dataf
% here we're reverse coding it in the variable "dv" from 0 - 4
% - by reverse coding, the smallest values in dv (i.e. the highest judgments of probability) 
% correspond to the smallest differences between last tone and predicted tone
% - minimum value is set to zero for model fitting purposes
dv = 5-dataf.resp_prob;

iv_list = {del_logf_pen,  del_logf_pred,  del_SD_pen,  del_SD_pred};
iv_labels2 = {'logf pen', 'logf pred', 'SD pen', 'SD pred'};
iv_labels = {'|log f_{t+1} - log f_t|', '|log f_{t+1} - E[log f_{t+1}]|', ...
             '|log f_{t+1} - log f_t| / \sigma_\epsilon', '|log f_{t+1} - E[log f_{t+1}]| / \sigma_\epsilon'};
% iv_labels2 = {'|log f_{t+1} - log f_t|', '|log f_{t+1} - E[log f_{t+1}|', ...
%              '|log f_{t+1} - log f_t| / s_e', '|log f_{t+1} - E[log f_{t+1}]| / s_e'};

         
% h=figure;
% fs = 15;
% set(h,'DefaultAxesFontSize',fs);
% for i_iv = 1:4
%     iv = iv_list{i_iv};
% 
%     subplot(2,2,i_iv);
%     plot(iv, dv, 'b.');
%     ylim([.5 5.5])
%     xlabel(iv_labels{i_iv});
%     ylabel('prob rating')
% 
% end


for i_iv = 1:4
    iv = iv_list{i_iv};
    
    m = [dv' iv'];
    in.DATA = m;
    out = MATLAB_Ordered_Probit_Estimate(in);

    logL(i_sub, i_iv) = out.Likelihood.LLV;

    % recover model probabilities
    x = iv*out.Beta;

    probit_x{i_iv, i_sub} = x;
    probit_beta{i_iv}(i_sub,:) = out.Beta;
    probit_cutpoints{i_iv,i_sub} = out.Cut_Points;

    % get the model-prediected final tone probability rating for each trial
    for i_trial = 1:length(x)
        % pm scales from 0 to 4
        % 0 --> final tone was close to prediction
        % 4 --> final tone was far from prediction
        pm = sum( x(i_trial) > out.Cut_Points );
        
        % now reverse-score pm so higher ratings -> higher tf probability
        prob_model{i_iv, i_sub}(i_trial) = 5 - pm;
    end
    
    
    
% % %     [beta mu l] = fit_ordered_probit_MLE(iv', dv', 4);
% % %     logL(i_sub, i_iv) = l;
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
% % %     logL(i_sub, i_iv) = l;
% % %     x = iv;
% % % 
% % %     probit_x{i_iv} = x;
% % %     probit_sd{i_iv} = sd;
% % %     probit_cutpoints{i_iv} = mu;
% % %     
% % %     for i_trial = 1:length(x)
% % %         prob_model{i_iv}(i_trial) = 1 + sum( x(i_trial) > mu );
% % %     end
    
    
    % make a heat map matrix for model prediction vs actual data
    for i_p = 1:5
        f = dataf.resp_prob == i_p;

        for i_pm = 1:5
            g = prob_model{i_iv, i_sub} == i_pm;
            prob_matrix{i_iv}(i_p, i_pm, i_sub) = sum(f & g) / sum(f);
            prob_matrix2{i_iv}(i_pm, i_p, i_sub) = sum(f & g) / length(f);
        end
    end
    
    [r_pmp(i_sub, i_iv) p_pmp(i_sub, i_iv)] = corr(prob_model{i_iv, i_sub}', dataf.resp_prob', 'type', 'Spearman');
    
    
%     for k = 1:5
%         p_b(k) = sum(dv==k)/length(dv);
%         p_m{i_iv}(k) = sum(prob_model{i_iv}==k)/length(prob_model{i_iv});
%     end
end


% get tone distances vs actual and modeled probability ratings
resp_probs{i_sub} = dataf.resp_prob;
tone_dists{i_sub} = [abs_pen; abs_pred; sd_pen; sd_pred];

for i_p = 1:5    
    f  = resp_probs{i_sub} == i_p;
    fm = prob_model{i_iv, i_sub} == i_p;
    for i_iv = 1:4
%         tone_dists_by_resp_prob_mean{i_sub}(i_iv, i_p) = mean(tone_dists{i_sub}(i_iv, f));
%         tone_dists_by_resp_prob_sem{i_sub}(i_iv, i_p) = std(tone_dists{i_sub}(i_iv, f)) / sqrt(length(tone_dists{i_sub}(i_iv, f)));

        td = Scale(tone_dists{i_sub}(i_iv, :));

        % tone dist by true prob rating
        tone_dists_by_resp_prob_mean{i_sub}(i_iv, i_p) = mean(td(f));
        tone_dists_by_resp_prob_sem{i_sub}(i_iv, i_p) = std(td(f)) / sqrt(length(td(f)));

        % tone dist by modeled prob rating
        tone_dists_by_resp_prob_mean_model{i_sub}(i_iv, i_p) = mean(td(fm));
        tone_dists_by_resp_prob_sem_model{i_sub}(i_iv, i_p) = std(td(fm)) / sqrt(length(td(fm)));

    end
end



% figure;
% subplot(211);
% b = [p_b' p_m{1}' p_m{2}' p_m{3}' p_m{4}'];
% bar(b);
% 
% subplot(212);
% bdiff = [p_b'-p_m{1}' p_b'-p_m{2}' p_b'-p_m{3}' p_b'-p_m{4}'];
% bar(bdiff);

K = 5; % ??
N = length(dv);

BIC = -2*logL(i_sub,:) + K.*log(N);
AIC = -2*logL(i_sub,:) + 2.*K.*N./(N-K-1);

BICw(i_sub,:) = exp( -.5*(BIC-min(BIC)) ) / sum( exp( -.5*(BIC-min(BIC)) ));
AICw(i_sub,:) = exp( -.5*(AIC-min(AIC)) ) / sum( exp( -.5*(AIC-min(AIC)) ));

% figure;
% fs = 15;
% 
% bar(AICw);
% set(gca, 'XTickLabel', iv_labels2);
% ylabel('Akaike weight', 'FontSize', fs)
% set(gca, 'FontSize', fs);

end


end   % end loop over subjects


%% post-subject loop analysis

%% adjust conf_br
% make it so that any (i_b, i_r) slot in conf_br that has a NaN for any
% subject has NaN for all subjects...

nanflag = zeros(5);
for i_b = 1:5
    for i_r = 1:5
        for i_s = 1:length(subs)
            if isnan(conf_br(i_r, i_b, i_s))
                nanflag(i_r, i_b) = 1;
            end
        end
    end
end


for i_s = 1:length(subs)
    xx = conf_br(:,:,i_s);
    xx(nanflag==1) = NaN;
    conf_br_nan(:,:,i_s) = xx;
end


%%

mean([r_abs_pen' r_abs_pred' r_sd_pen' r_sd_pred'])

mean(r_abs_pen_b)
mean(r_abs_pred_b)
mean(r_sd_pen_b)
mean(r_sd_pred_b)


%% beta discrimination vs final tone prediction

type = 'Pearson';
[r_beta_abs_pen  p_beta_abs_pen ] = corr(r_beta', r_abs_pen',  'type', type)
[r_beta_abs_pred p_beta_abs_pred] = corr(r_beta', r_abs_pred', 'type', type)
[r_beta_sd_pen   p_beta_sd_pen  ] = corr(r_beta', r_sd_pen',   'type', type)
[r_beta_sd_pred  p_beta_sd_pred ] = corr(r_beta', r_sd_pred',  'type', type)


%% plot data

figure;
fs = 15;

imagesc(mean(beta_matrix,3))
axis square
colormap hot
colorbar

xlabel('response beta', 'FontSize', fs)
ylabel('stimulus beta', 'FontSize', fs)
title(['mean p( response beta | stim beta ), rho = ' num2str( mean(r_beta) )], 'FontSize', fs)

set(gca,'XTick',1:5)
set(gca,'YTick',1:5)
set(gca,'XTickLabel',{'0' '0.5' '1' '1.5' '2'})
set(gca,'YTickLabel',{'0' '0.5' '1' '1.5' '2'})
set(gca,'FontSize',fs)

% caxis([0 .71])




figd; hold on;

for i_sub = 1:nSubs
    beta_d_s = beta_d(:,:,i_sub);
    beta_d_s(abs(beta_d_s)==Inf) = NaN;
    dd = nanmean(beta_d_s,2);
    dcum(i_sub,:) = cumsum(dd');
end

errorbar([0 .5:.5:2], [0 mean(dcum)], [0 std(dcum)]./sqrt(nSubs), 'bo');
% plot([0 .5:.5:2], [0 dcum], 'bo');
pp = polyfit(.5:.5:2, mean(dcum), 1);
xx = [.5 2];
plot(xx, xx*pp(1)+pp(2), 'b-'); 
xlabel('beta')
ylabel('cumulative d''');
xlim([-.2 2.2])


%% confidence plots

% conf | accuracy as a function of beta
figd; hold on;
errorbar([0 .5:.5:2], mean(conf), std(conf)./sqrt(nSubs), 'ko-');
errorbar([0 .5:.5:2], mean(conf_corr0), std(conf_corr0)./sqrt(nSubs), 'rv-');
errorbar([0 .5:.5:2], mean(conf_corr1), std(conf_corr1)./sqrt(nSubs), 'g^-');
xlabel('stimulus beta')
ylabel('confidence in beta response')
legend('overall', 'incorrect', 'correct', 'Location', 'NorthWest')


% conf | error as a function of beta
figd; hold on;

cd_m = mean(conf_diff,3);
cd_s = std(conf_diff,0,3);

errorbar([0 .5:.5:2], cd_m(1,:), cd_s(1,:)./sqrt(nSubs), 'ro-');
errorbar([0 .5:.5:2], cd_m(2,:), cd_s(2,:)./sqrt(nSubs), 'ko-');
errorbar([0 .5:.5:2], cd_m(3,:), cd_s(3,:)./sqrt(nSubs), 'go-');

xlabel('stimulus beta')
ylabel('confidence in beta response')
legend('beta error >= 1', 'beta error = .5', 'beta error = 0', 'Location', 'NorthWest')


% % % % conf | error as a function of d'
% % % figd; hold on;
% % % 
% % % cd_m = mean(conf_diff,3);
% % % cd_s = std(conf_diff,0,3);
% % % 
% % % % errorbar([0 mean(dcum)], cd_m(1,:), cd_s(1,:)./sqrt(nSubs), 'ro-');
% % % % errorbar([0 mean(dcum)], cd_m(2,:), cd_s(2,:)./sqrt(nSubs), 'ko-');
% % % % errorbar([0 mean(dcum)], cd_m(3,:), cd_s(3,:)./sqrt(nSubs), 'go-');
% % % 
% % % errorbar(mean(dcum), cd_m(1,2:end), cd_s(1,2:end)./sqrt(nSubs), 'ro-');
% % % errorbar(mean(dcum), cd_m(2,2:end), cd_s(2,2:end)./sqrt(nSubs), 'ko-');
% % % errorbar(mean(dcum), cd_m(3,2:end), cd_s(3,2:end)./sqrt(nSubs), 'go-');
% % % 
% % % xlabel('cumulative d''')
% % % ylabel('confidence in beta response')
% % % legend('beta error >= 1', 'beta error = .5', 'beta error = 0', 'Location', 'NorthWest')


% conf | stim, response
figd; hold on;
cbr_m = mean(conf_br_nan,3);
cbr_s = std(conf_br_nan,0,3);

errorbar([0 .5:.5:2], cbr_m(1,:), cbr_s(1,:)./sqrt(nSubs), 'ko-');
errorbar([0 .5:.5:2], cbr_m(2,:), cbr_s(2,:)./sqrt(nSubs), 'mo-');
errorbar([0 .5:.5:2], cbr_m(3,:), cbr_s(3,:)./sqrt(nSubs), 'ro-');
errorbar([0 .5:.5:2], cbr_m(4,:), cbr_s(4,:)./sqrt(nSubs), 'go-');
errorbar([0 .5:.5:2], cbr_m(5,:), cbr_s(5,:)./sqrt(nSubs), 'bo-');

xlabel('stimulus beta')
ylabel('confidence in beta response')
legend('beta resp = 0', 'beta resp = .5', 'beta resp = 1', 'beta resp = 1.5', 'beta resp = 2', 'Location', 'NorthWest')



figd; hold on;
cbr_m2 = cbr_m';
cbr_s2 = cbr_s';

errorbar([0 .5:.5:2], cbr_m2(1,:), cbr_s2(1,:)./sqrt(nSubs), 'ko-');
errorbar([0 .5:.5:2], cbr_m2(2,:), cbr_s2(2,:)./sqrt(nSubs), 'mo-');
errorbar([0 .5:.5:2], cbr_m2(3,:), cbr_s2(3,:)./sqrt(nSubs), 'ro-');
errorbar([0 .5:.5:2], cbr_m2(4,:), cbr_s2(4,:)./sqrt(nSubs), 'go-');
errorbar([0 .5:.5:2], cbr_m2(5,:), cbr_s2(5,:)./sqrt(nSubs), 'bo-');

xlabel('response beta')
ylabel('confidence in beta response')
legend('beta stim = 0', 'beta stim = .5', 'beta stim = 1', 'beta stim = 1.5', 'beta stim = 2', 'Location', 'NorthWest')


%% confidence matrix

% overall conf matrix
figure;
fs = 15;

imagesc(mean(conf_matrix,3))
axis square
colormap hot
colorbar

xlabel('confidence in beta response', 'FontSize', fs)
ylabel('error in beta response', 'FontSize', fs)
title(['mean conf | error, rho = ' num2str( mean(r_conf) )], 'FontSize', fs)

set(gca,'XTick',1:5)
set(gca,'YTick',1:5)
set(gca,'XTickLabel',{'0' '1' '2' '3' '4'})
set(gca,'YTickLabel',{'>= 1' '0.5' '0'})
set(gca,'FontSize',fs)



% conf matrix for each level of beta
figure;
fs = 15;

i_p = 0;
for i_b = 1:5

    i_p = i_p + 1;
    
    subplot(3,2,i_p);
    imagesc(mean(conf_matrix_b{i_b},3))
%     axis square
    colormap hot
    colorbar

    xlabel('confidence in beta response', 'FontSize', fs)
    ylabel('error in beta response', 'FontSize', fs)

    title({'mean conf | error', ...
        ['true beta = ' num2str(betas(i_b)) ', rho = ' num2str( mean(r_conf_b{i_b}) )]}, 'FontSize', fs)

    set(gca,'XTick',1:5)
    set(gca,'YTick',1:5)
    set(gca,'XTickLabel',{'0' '1' '2' '3' '4'})
    set(gca,'YTickLabel',{'>= 1' '0.5' '0'})
    set(gca,'FontSize',fs)

end





% conf matrix for each level of reported beta
figure;
fs = 15;

i_p = 0;
for i_r = 1:5

    i_p = i_p + 1;
    
    subplot(3,2,i_p);
    imagesc(mean(conf_matrix_r{i_r},3))
%     axis square
    colormap hot
    colorbar

    xlabel('confidence in beta response', 'FontSize', fs)
    ylabel('error in beta response', 'FontSize', fs)

    title({'mean conf | error', ...
        ['reported beta = ' num2str(betas(i_r)) ', rho = ' num2str( mean(r_conf_r{i_r}) )]}, 'FontSize', fs)

    set(gca,'XTick',1:5)
    set(gca,'YTick',1:5)
    set(gca,'XTickLabel',{'0' '1' '2' '3' '4'})
    set(gca,'YTickLabel',{'>= 1' '0.5' '0'})
    set(gca,'FontSize',fs)

end


%% final tone probability rating

% probability rating vs beta discrimination
figd;

for i_p = 1:4
    subplot(2,2,i_p); hold on;
    xlabel('rho(stim beta, resp beta)');
    ylabel('r(prob rating, tone dist)');
    
    if i_p == 1
        v = r_abs_pen;
    elseif i_p == 2
        v = r_abs_pred;
    elseif i_p == 3
        v = r_sd_pen;
    elseif i_p == 4
        v = r_sd_pred;
    end
        
    plot(r_beta, v, 'bo');
    pp = polyfit(r_beta, v, 1);
    xlim([.3 .8])
    xx = xlim;
    yy = ylim;
    plot(xx, xx*pp(1)+pp(2), 'b-'); 
    
    [r p] = corr(r_beta', v',  'type', 'Pearson');
    title({ iv_labels{i_p}, ...
        ['r = ' num2str(r,2) ', p = ' num2str(p,2)]});
end

        
        
        
% probability rating vs beta level
iv_labels = {'|log f_{t+1} - log f_t|', '|log f_{t+1} - E[log f_{t+1}]|', ...
             '|log f_{t+1} - log f_t| / \sigma_\epsilon', '|log f_{t+1} - E[log f_{t+1}]| / \sigma_\epsilon'};

figd; hold on;

r_m = mean(r_abs_pen_b);
r_s = std(r_abs_pen_b);
errorbar([0 .5:.5:2], r_m(1,:), r_s(1,:)./sqrt(nSubs), 'bo-');

r_m = mean(r_abs_pred_b);
r_s = std(r_abs_pred_b);
errorbar([0 .5:.5:2], r_m(1,:), r_s(1,:)./sqrt(nSubs), 'ro-');

r_m = mean(r_sd_pen_b);
r_s = std(r_sd_pen_b);
errorbar([0 .5:.5:2], r_m(1,:), r_s(1,:)./sqrt(nSubs), 'ko-');

r_m = mean(r_sd_pred_b);
r_s = std(r_sd_pred_b);
errorbar([0 .5:.5:2], r_m(1,:), r_s(1,:)./sqrt(nSubs), 'go-');

xlabel('stimulus beta')
ylabel('Pearson''s r')
title('correlation for final tone prob rating and final tone dist from reference')
% legend('beta resp = 0', 'beta resp = .5', 'beta resp = 1', 'beta resp = 1.5', 'beta resp = 2', 'Location', 'NorthWest')
legend(iv_labels{1}, iv_labels{2}, iv_labels{3}, iv_labels{4}, 'Location', 'NorthWest')
set(gca,'XTick', 0:.5:2);


%% AIC analysis

if model_selection

    
% % plot of Akaike weights for each subject
% figure;

    
% plot of actual vs model-predicted final tone probability ratings
for i_sub = 1:length(subs)


figd; 
for i_iv = 1:4
    subplot(3,4,i_iv); hold on;
    
    tone_dist_mean = tone_dists_by_resp_prob_mean{i_sub}(i_iv,:);
    tone_dist_sem  = tone_dists_by_resp_prob_sem{i_sub}(i_iv,:);
    
    tone_dist_mean_model = tone_dists_by_resp_prob_mean_model{i_sub}(i_iv,:);
    tone_dist_sem_model  = tone_dists_by_resp_prob_sem_model{i_sub}(i_iv,:);    
    
    % get the marginal probability for the subject giving each probability rating
    pr_prob   = sum( prob_matrix2{i_iv}(:,:,i_sub) );
    pr_prob_m = sum( prob_matrix2{i_iv}(:,:,i_sub), 2);
    
    % plot the data
    offset = .2; width = .35; lw = 3;
    for i_p = 1:5
        % behavioral
        h = bar(i_p - offset, tone_dist_mean(i_p), width);
        set(h, 'FaceColor', (1-pr_prob(i_p))*ones(1,3));
        set(h, 'EdgeColor', [1 0 0], 'LineWidth', lw);
        
        % modeled
        h = bar(i_p + offset, tone_dist_mean_model(i_p), width);
        set(h, 'FaceColor', (1-pr_prob_m(i_p))*ones(1,3));
        set(h, 'EdgeColor', [0 0 1], 'LineWidth', lw);        
    end
    ylim([0 1]);
    
    % h = bar([tone_dist_mean' tone_dist_mean_model']);

    set(gca, 'XTick', 1:5);
    set(gca, 'XTickLabel', {'1', '2', '3', '4', '5'});
    
    if i_iv == 1
        legend('behav', 'modeled', 'Location', 'NorthEast');
        xlabel('prob rating')
        ylabel('normalized tone dist')
    end
    
    title(iv_labels{i_iv});

    
    
% % %     h = bar([tone_dist_mean' tone_dist_mean_model']);
% % % %     errorbar(1:5, tone_dist_mean, tone_dist_sem, 'k.')
% % %     if i_iv == 1
% % %         xlabel('prob rating')
% % %         ylabel('tone dist')
% % %     end
% % %     set(gca, 'XTick', 1:5);
% % %     set(gca, 'XTickLabel', {'1', '2', '3', '4', '5'});

    
    h=subplot(3,4,i_iv+4);
    imagesc(prob_matrix2{i_iv}(:,:,i_sub)); 
    colormap(hot);

    caxis([0 .26]);
    b = mycolorbar;
    set(b, 'YTick', [0 .1 .2]);
    
    axis square
    set(gca, 'XTick', 1:5);
    set(gca, 'XTickLabel', {'1', '2', '3', '4', '5'});
    set(gca, 'YTick', 1:5);
    set(gca, 'YTickLabel', {'1', '2', '3', '4', '5'});
    
    
%     if i_iv==4
%         b=colorbar;
%         caxis([0 .26]);
%         pos = get(h, 'Position');
% %         keyboard
%         set(b,'Position', pos.*[.95 .12 .05 .8]);
%     end
    
    if i_iv == 1
        ylabel('prob rating (model)');
        xlabel('prob rating (behav)');
    elseif i_iv == 2
        title('joint probability distributions...');
    elseif i_iv == 3
        title('...for behav & model prob rating');
    end
    
end % end i_iv loop

subplot(3,4,9:12);

% modelStats1 = [r_abs_pen(i_sub),  r_pmp(i_sub,1), AICw(i_sub,1)];
% modelStats2 = [r_abs_pred(i_sub), r_pmp(i_sub,2), AICw(i_sub,2)];
% modelStats3 = [r_sd_pen(i_sub),   r_pmp(i_sub,3), AICw(i_sub,3)];
% modelStats4 = [r_sd_pred(i_sub),  r_pmp(i_sub,4), AICw(i_sub,4)];

r_perf  = [r_abs_pen(i_sub), r_abs_pred(i_sub), r_sd_pen(i_sub), r_sd_pred(i_sub)];
r_model = r_pmp(i_sub,:);
Akaike_w = AICw(i_sub, :);

bar([r_perf' r_model' Akaike_w'])
title(['Subject ' subs{i_sub}]);

set(gca,'XTickLabel',{})
legend('rho(tone dist/behav prob rating)', 'rho(modeled/behav prob rating)', 'Akaike weight')

% mtit(['Subject ' subs{i_sub} ', p(model rating & true rating)'], 'xoff', -.1);

end % end i_sub loop

end
    
%     for i_trial = 1:length(x)
%         prob_model{i_iv}(i_sub, i_trial) = 1 + sum( x(i_trial) > out.Cut_Points );
%     end