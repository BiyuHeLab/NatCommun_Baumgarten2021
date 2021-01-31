function  NASTD_MEG_Behav_FTPLStat_Subs...
    (subs, tonedur_text, ...
    paths_NASTD_MEG)
%Compute single subject statistics for FTPL ratings:
%Computes single-subject ANOVA to get F-value quantifying the interaction
%effect of predp34 and p34 final tone on FTPL rating.

%Performs 2 ANOVAs for each subject seperately, in order to get F-values for each subject
%1. 3way ANOVA, factors: math. predp34 pitch, p34 pitch, tone duration - output = Fsub.mat (i.e., across tone durations)
%2. 2way ANOVA, factors: math. predp34 pitch, p34 pitch - output = FsubDur.mat (i.e., within tone durations)

%% 0) Specify paths and analysis options
path_outputdata = ([paths_NASTD_MEG.Current_outputdata 'Statistics/']);
mkdir(path_outputdata);

%% 1) load behavioral data and select specific trials
%Bring behavioral data into format for ANOVA
[subPred, subDurPred] = NASTD_MEG_Behav_FTPL_PrepData4ANOVA(subs, tonedur_text, paths_NASTD_MEG);

%Delete NaN entries (no response/missed trials)
for i_sub = 1:length(subs)
    NAN_trials = [];
    NAN_trials = find(isnan(subPred{i_sub}.PredRating));
    if ~isempty(NAN_trials)
        disp([num2str(length(NAN_trials)) ' NaN trials found for sub ' num2str(i_sub)])
        for i_NANtrials = NAN_trials
            subPred{i_sub}.PredRating(i_NANtrials) = [];
            subPred{i_sub}.tonedur(i_NANtrials) = [];
            subPred{i_sub}.predp34(i_NANtrials) = [];
            subPred{i_sub}.p34(i_NANtrials) = [];
        end
    end
end

%% 2) Compute Statistics for FTPL ratings:
% 2.1)  3way ANOVA (dependent var: FTPL rating;
%factors: tone duration, math. predp34 FTP, presented FTP
varnames = {'tonedur'; 'predp34'; 'p34'};
for i_sub = 1:length(subs)
    
    [~, tbl{i_sub}] = anovan(subPred{i_sub}.PredRating, {subPred{i_sub}.tonedur subPred{i_sub}.predp34 subPred{i_sub}.p34},...
        'model', 'interaction', 'varnames', varnames, 'display', 'off');
    
    FTPLstat.SingleSub.ANOVA3.FP_Main_TD(i_sub,:) = tbl{i_sub}(2,6:7); %Main effect for tone duration
    FTPLstat.SingleSub.ANOVA3.FP_Main_p34(i_sub,:) = tbl{i_sub}(4,6:7); %Main effect for p34
    FTPLstat.SingleSub.ANOVA3.FP_Main_predp34(i_sub,:) = tbl{i_sub}(3,6:7);%Main effect for p*34    
    FTPLstat.SingleSub.ANOVA3.FP_Inter_p34predp34(i_sub,:) = tbl{i_sub}(7,6:7);%Predp34*p34 pitch interaction effect

end

% 2.2)  2way ANOVA - 1 ANOVA per tone duration
%(dependent var: final tone likelihood rating;
%factors: math. predp34 FTP, p34 FTP,
tone_durs = [0.15 0.3 0.6];
varnames2 = {'predp34'; 'p34'};

for i_tonedur = 1:length(tone_durs)
    for i_sub = 1:length(subs)
        
        [~, subDurTbl{i_sub,i_tonedur}] = ...
            anovan(subDurPred{i_sub,i_tonedur}.PredRating, {subDurPred{i_sub,i_tonedur}.predp34 subDurPred{i_sub,i_tonedur}.p34},...
            'model', 'interaction', 'varnames', varnames2, 'display', 'off');
        
        FsubDur(i_sub,i_tonedur) = subDurTbl{i_sub,i_tonedur}(4,6);
        %Subject- & toneduration-wise F value for math.predp34*presented final tone pitch interaction effect
        
        FTPLstat.SingleSub.ANOVA2.FP_Main_p34{i_sub,i_tonedur} = subDurTbl{i_sub, i_tonedur}(3,6:7);%Main effect for p34
        FTPLstat.SingleSub.ANOVA2.FP_Main_predp34{i_sub,i_tonedur} = subDurTbl{i_sub, i_tonedur}(2,6:7);%Main effect for p*34    
        FTPLstat.SingleSub.ANOVA2.FP_Inter_p34predp34{i_sub,i_tonedur} = subDurTbl{i_sub, i_tonedur}(4,6:7);%Predp34*p34 pitch interaction effect
        
    end
end

%% 3) Save data
save([path_outputdata 'FTPLrating_Stats_Subs.mat'], 'FTPLstat');
