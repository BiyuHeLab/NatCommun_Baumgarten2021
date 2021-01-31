function  NASTD_MEG_Behav_FTPLStat_Group...
    (subs, tonedur_text, ...
    paths_NASTD_MEG)

%Compute group level behavioral Final Tone Pitch Likelihood (FTPL) rating 
%effect across and for each tone duration

%1) Computes group level repeated measures ANOVAs with DV FTPL ratings
%IVs: p34,p*34, Tone Durations
%IVs: p34,p*34 (per Tone Duration)

%2) Computes group level t-proxy_matrixs bewteen extreme math. predp34 - presented
%FTP (i.e., math.predp34 FTP = low + presented FTP = low vs. math.predp34
%FTP = high + presened FTP = low)

%Optional: Restrict trial selection to first or last half of the experiment

%% 0) Specify paths and analysis options
path_outputdata = ([paths_NASTD_MEG.Current_outputdata 'Statistics/']); %Output: Fvals
mkdir(path_outputdata);

%% 1) load behavioral data and select specific trials
%Bring behavioral data into format for ANOVA
[subPred, subDurPred] = NASTD_MEG_Behav_FTPL_PrepData4ANOVA(subs, tonedur_text, paths_NASTD_MEG);

% %Optional: Restrict trials to first/last half of recording
% %Across TD
% for i_sub = 1:length(subs)
%     
%     subPred_SelTrials{i_sub}.PredRating = [];
%     subPred_SelTrials{i_sub}.tonedur_numerical = [];
%     subPred_SelTrials{i_sub}.predp34_numerical = [];
%     subPred_SelTrials{i_sub}.p34_numerical = [];
%     subPred_SelTrials{i_sub}.beta_numerical = [];
%     
%     proxy_matrix = [subPred{i_sub}.PredRating, ...
%         subPred{i_sub}.tonedur_numerical, ...
%         subPred{i_sub}.predp34_numerical, ...
%         subPred{i_sub}.p34_numerical, ...
%         subPred{i_sub}.beta_numerical];
%     
%     allselTrials{i_sub} = [];
%         
%     for i_tonedur = 1:length(unique(subPred{i_sub}.tonedur))
%         for i_predp34 = 1:3
%             for i_p34 = 1:6
%                 for i_beta = 1:3
%                     
%                     %Trial selection based on TD, predp34, p34, and beta
%                     selected_rows = find(proxy_matrix(:,2) == i_tonedur & ...
%                         proxy_matrix(:,3) == i_predp34 & ...
%                         proxy_matrix(:,4) == i_p34 & ...
%                         proxy_matrix(:,5) == i_beta);
%                     
%                     selTrials = selected_rows(1:(round(length(selected_rows)/2))); %first half of experiment
%                     % selTrials = selected_rows(round(length(selected_rows)/2)+1:length(selected_rows)); %second half of experiment
%                     
%                     subPred_SelTrials{i_sub}.PredRating =  ...
%                         [subPred_SelTrials{i_sub}.PredRating; proxy_matrix(selTrials,1)];
%                     subPred_SelTrials{i_sub}.tonedur_numerical =  ...
%                         [subPred_SelTrials{i_sub}.tonedur_numerical; proxy_matrix(selTrials,2)];
%                     subPred_SelTrials{i_sub}.predp34_numerical = ...
%                         [subPred_SelTrials{i_sub}.predp34_numerical; proxy_matrix(selTrials,3)] ;
%                     subPred_SelTrials{i_sub}.p34_numerical = ...
%                         [subPred_SelTrials{i_sub}.p34_numerical; proxy_matrix(selTrials,4)];
%                     subPred_SelTrials{i_sub}.beta_numerical = ...
%                         [subPred_SelTrials{i_sub}.beta_numerical; proxy_matrix(selTrials,5)];
%                     allselTrials{i_sub} = [allselTrials{i_sub}; selTrials];
%                 end
%             end
%         end
%     end
% end
% 
% %Per TD
% for i_sub = 1:length(subs)
%     for i_tonedur = 1:length(unique(subPred{i_sub}.tonedur))
%         
%         subDurPred_SelTrials{i_sub,i_tonedur}.PredRating = [];
%         subDurPred_SelTrials{i_sub,i_tonedur}.predp34_numerical = [];
%         subDurPred_SelTrials{i_sub,i_tonedur}.p34_numerical = [];
%         subDurPred_SelTrials{i_sub,i_tonedur}.beta_numerical = [];
%         
%         allselTrials_perTD{i_sub,i_tonedur} = [];
%         
%         proxy_matrix_perTD = [subDurPred{i_sub,i_tonedur}.PredRating, ...
%             subDurPred{i_sub,i_tonedur}.predp34_numerical, ...
%             subDurPred{i_sub,i_tonedur}.p34_numerical, ...
%             subDurPred{i_sub,i_tonedur}.beta_numerical];
%         
%         for i_predp34 = 1:3
%             for i_p34 = 1:6
%                 for i_beta = 1:3
%                     
%                     selected_rows = find(proxy_matrix_perTD(:,2) == i_predp34 & ...
%                         proxy_matrix_perTD(:,3) == i_p34 & ...
%                         proxy_matrix_perTD(:,4) == i_beta);
%                     
%                     selTrials = selected_rows(1:(round(length(selected_rows)/2))); %first half of experiment
%                     %  selTrials =  selected_rows(round(length(selected_rows)/2)+1:length(selected_rows)); %second half of experiment
%                     
%                     subDurPred_SelTrials{i_sub,i_tonedur}.PredRating =  ...
%                         [subDurPred_SelTrials{i_sub,i_tonedur}.PredRating; proxy_matrix_perTD(selTrials,1)];
%                     subDurPred_SelTrials{i_sub,i_tonedur}.predp34_numerical = ...
%                         [subDurPred_SelTrials{i_sub,i_tonedur}.predp34_numerical; proxy_matrix_perTD(selTrials,2)] ;
%                     subDurPred_SelTrials{i_sub,i_tonedur}.p34_numerical = ...
%                         [subDurPred_SelTrials{i_sub,i_tonedur}.p34_numerical; proxy_matrix_perTD(selTrials,3)];
%                     subDurPred_SelTrials{i_sub,i_tonedur}.beta_numerical = ...
%                         [subDurPred_SelTrials{i_sub,i_tonedur}.beta_numerical; proxy_matrix_perTD(selTrials,4)];
%                     
%                     allselTrials_perTD{i_sub, i_tonedur} = [allselTrials_perTD{i_sub, i_tonedur}; selTrials];
%                 end
%             end
%         end
%     end
% end
% 
% subPred = subPred_SelTrials;
% subDurPred = subDurPred_SelTrials;

%% 2) Group level FTPL rating comparisons via ANOVA:
%2.1 across Tonedur: 3-way RM ANOVA (predp34, p34, ToneDur)
%Average FTPL ratings across trials with identical  predp34, p34, and subject
FTPLrating_AvgTrials = [];
predp34_num_AvgTrials = [];
p34_num_AvgTrials = [];
tonedur_num_AvgTrials = [];
subjects_AvgTrials = [];

for i_sub = 1:length(subs) %Concatenate all trials for all subs per tone condition and fill proxies
    proxy_matrix = [subPred{i_sub}.PredRating, ...
        subPred{i_sub}.tonedur_numerical, ...
        subPred{i_sub}.predp34_numerical, ...
        subPred{i_sub}.p34_numerical,...
        subPred{i_sub}.beta_numerical];
    
    for i_tonedur = 1:3
        for i_predp34 = 1:3
            for i_p34 = 1:6
                
                selected_rows = find(proxy_matrix(:,2) == i_tonedur & ...
                    proxy_matrix(:,3) == i_predp34 & ...
                    proxy_matrix(:,4) == i_p34);
                
                FTPLrating_AvgTrials = ...
                    [FTPLrating_AvgTrials; nanmean(proxy_matrix(selected_rows,1))];
                tonedur_num_AvgTrials = ...
                    [tonedur_num_AvgTrials; repelem(i_tonedur,1)'];
                predp34_num_AvgTrials = ...
                    [predp34_num_AvgTrials; repelem(i_predp34,1)'];
                p34_num_AvgTrials = ...
                    [p34_num_AvgTrials; repelem(i_p34,1)'];
                subjects_AvgTrials = ...
                    [subjects_AvgTrials; repelem(i_sub,1)'];
                
            end
        end
    end
    proxy_matrix = [];
end

%Delete NaN entries (no response/missed trials)
NAN_trials = [];
NAN_trials = find(isnan(FTPLrating_AvgTrials));
if ~isempty(NAN_trials)
    for i_NANtrials = NAN_trials
        FTPLrating_AvgTrials(i_NANtrials) = [];
        tonedur_num_AvgTrials(i_NANtrials) = [];
        predp34_num_AvgTrials(i_NANtrials) = [];
        p34_num_AvgTrials(i_NANtrials) = [];
        subjects_AvgTrials(i_NANtrials) = [];
    end
end

[~,table_ANOVA3,~] = anovan(FTPLrating_AvgTrials,...
    {tonedur_num_AvgTrials,...
    predp34_num_AvgTrials,...
    p34_num_AvgTrials,...
    subjects_AvgTrials},...
    'model', 'full','random',4,'varnames',{'ToneDur','p*34','p34','subjects'},'display','off');

%2.2 per ToneDuration: 2-way RM ANOVA (predp34, p34)
tone_durs = [0.15 0.3 0.6];
%Concatenate trials across subjects per tone duration condition and compute
for i_tonedur = 1:length(tone_durs) %Separate RM ANOVA for each tone duration

    %Average FTPL ratings across trials with identical  predp34, p34, and subject
    FTPLrating_AvgTrials_perTD{i_tonedur} = [];
    predp34_num_AvgTrials_perTD{i_tonedur} = [];
    p34_num_AvgTrials_perTD{i_tonedur} = [];
    subjects_AvgTrials_perTD{i_tonedur} = [];

    for i_sub = 1:length(subs) %Concatenate all trials for all subs per tone condition and fill proxies    
        proxy_matrix = [subDurPred{i_sub,i_tonedur}.PredRating, ...
            subDurPred{i_sub,i_tonedur}.predp34_numerical, ...
            subDurPred{i_sub,i_tonedur}.p34_numerical];

        for i_predp34 = 1:3
            for i_p34 = 1:6
                selected_rows = ...
                    find(proxy_matrix(:,2) == i_predp34 & proxy_matrix(:,3) == i_p34);            
                FTPLrating_AvgTrials_perTD{i_tonedur} = ...
                    [FTPLrating_AvgTrials_perTD{i_tonedur}; nanmean(proxy_matrix(selected_rows,1))];
                predp34_num_AvgTrials_perTD{i_tonedur} = ...
                    [predp34_num_AvgTrials_perTD{i_tonedur}; nanmean(proxy_matrix(selected_rows,2))];
                p34_num_AvgTrials_perTD{i_tonedur} = ...
                    [p34_num_AvgTrials_perTD{i_tonedur}; nanmean(proxy_matrix(selected_rows,3))];
                subjects_AvgTrials_perTD{i_tonedur} = ...
                    [subjects_AvgTrials_perTD{i_tonedur}; repelem(i_sub,1)'];
            end
        end
        proxy_matrix = [];
    end    
    
    %Delete NaN entries (no response/missed trials)
    NAN_trials = [];
    NAN_trials = find(isnan(FTPLrating_AvgTrials_perTD{i_tonedur}));
    if ~isempty(NAN_trials)
        for i_NANtrials = NAN_trials
            FTPLrating_AvgTrials_perTD{i_tonedur}(i_NANtrials) = [];
            predp34_num_AvgTrials_perTD{i_tonedur}(i_NANtrials) = [];
            p34_num_AvgTrials_perTD{i_tonedur}(i_NANtrials) = [];
            subjects_AvgTrials_perTD{i_tonedur}(i_NANtrials) = [];
        end
    end     
    
    %Repeated Measures Two-way Analysis of Variance proxy_matrix
    %dependent var: final tone likelihood rating;
    %factors/iv's: math. predp34 FTP, presented FTP
    [~,table_ANOVA2{i_tonedur},~] = ...
        anovan(FTPLrating_AvgTrials_perTD{i_tonedur},...
        {p34_num_AvgTrials_perTD{i_tonedur},...
        predp34_num_AvgTrials_perTD{i_tonedur},...
        subjects_AvgTrials_perTD{i_tonedur}},...
        'model', 'full','random',3,'varnames',{'p34','p*34','subjects'},'display','off');
    
end

%% 3) Group level comparison via t-proxy_matrix
%Compare maxima (low & high likelihood ratings) across subjects
for i_sub = 1:length(subs)
    for i_tonedur = 1:length(tone_durs) %no all conditions
        %Select trials with extrema in math. predp34 (low vs high) and presented FTPs (-3 vs 3)
        f_expLow = strcmp(subDurPred{i_sub,i_tonedur}.predp34,'low');
        f_expHigh = strcmp(subDurPred{i_sub,i_tonedur}.predp34,'high');
        f_actLow = strcmp(subDurPred{i_sub,i_tonedur}.p34,'-3');
        f_actHigh = strcmp(subDurPred{i_sub,i_tonedur}.p34,'3');
        
        %Compute average rating for extrema trials for each subject for each tone duration
        mu_actLow_expLow_perTonDur(i_sub,i_tonedur) = ...
            nanmean(subDurPred{i_sub,i_tonedur}.PredRating(f_actLow & f_expLow));
        mu_actLow_expHigh_perTonDur(i_sub,i_tonedur) = ...
            nanmean(subDurPred{i_sub,i_tonedur}.PredRating(f_actLow & f_expHigh));
        mu_actHigh_expLow_perTonDur(i_sub,i_tonedur) = ...
            nanmean(subDurPred{i_sub,i_tonedur}.PredRating(f_actHigh & f_expLow));
        mu_actHigh_expHigh_perTonDur(i_sub,i_tonedur) = ...
            nanmean(subDurPred{i_sub,i_tonedur}.PredRating(f_actHigh & f_expHigh));
    end
    
    %Average over tone durations and place average in 4th column
    mu_actLow_expLow_perTonDur(i_sub,4) = ...
        nanmean(mu_actLow_expLow_perTonDur(i_sub,:));
    mu_actLow_expHigh_perTonDur(i_sub,4) = ...
        nanmean(mu_actLow_expHigh_perTonDur(i_sub,:));
    mu_actHigh_expLow_perTonDur(i_sub,4) = ...
        nanmean(mu_actHigh_expLow_perTonDur(i_sub,:));
    mu_actHigh_expHigh_perTonDur(i_sub,4)= ...
        nanmean(mu_actHigh_expHigh_perTonDur(i_sub,:));
end

for i_tonedur = 1:size(mu_actHigh_expHigh_perTonDur,2)
    [h_low, p_low(i_tonedur),~, stat_low(i_tonedur)] = ...
        ttest2(mu_actLow_expLow_perTonDur(:,i_tonedur), mu_actLow_expHigh_perTonDur(:,i_tonedur));
    diff_low(i_tonedur) = ...
        mean(mu_actLow_expLow_perTonDur(:,i_tonedur)) - ...
        mean(mu_actLow_expHigh_perTonDur(:,i_tonedur));
    CohensD_difflow(i_tonedur) = ...
        mean(mu_actLow_expLow_perTonDur(:,i_tonedur) - ...
        mu_actLow_expHigh_perTonDur(:,i_tonedur)) / ...
        std(mu_actLow_expLow_perTonDur(:,i_tonedur) - ...
        mu_actLow_expHigh_perTonDur(:,i_tonedur));
    
    [h_high,p_high(i_tonedur),~, stat_high(i_tonedur)] = ...
        ttest2(mu_actHigh_expLow_perTonDur(:,i_tonedur), ...
        mu_actHigh_expHigh_perTonDur(:,i_tonedur));
    diff_high(i_tonedur) = ...
        mean(mu_actHigh_expLow_perTonDur(:,i_tonedur)) - ...
        mean(mu_actHigh_expHigh_perTonDur(:,i_tonedur));
    CohensD_diffhigh(i_tonedur) = ...
        mean(mu_actHigh_expLow_perTonDur(:,i_tonedur) - ...
        mu_actHigh_expHigh_perTonDur(:,i_tonedur)) / ...
        std(mu_actHigh_expLow_perTonDur(:,i_tonedur) - ...
        mu_actHigh_expHigh_perTonDur(:,i_tonedur));
end

%% Optional: Bring FTPL data into SPSS-conform formet and save as .xls
% CombMatrix = ...
%     [FTPLrating_AvgTrials, ...
%     tonedur_num_AvgTrials, ...
%     predp34_num_AvgTrials, ...
%     p34_num_AvgTrials, ...
%     subjects_AvgTrials];
% 
% counter_column = 1;
% 
% for i_tonedur = 1:3
%     for i_predp34 = 1:3
%         for i_p34 = 1:6
%             ind = find(CombMatrix(:,2) == i_tonedur & CombMatrix(:,3) == i_predp34 & CombMatrix(:,4) == i_p34);
%             SPSSdata(:, counter_column) = CombMatrix(ind,1);
%             counter_column = counter_column + 1;
%         end
%     end
% end
% 
% %Save output matrix as Excel file, which can be read into SPSS
% cd(paths_NASTD_MEG.Analysis.Behavior)
% xlswrite('SPSSdata.xlsx',SPSSdata)

end
