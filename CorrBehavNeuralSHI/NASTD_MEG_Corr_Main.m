%% Baumgarten et al. Nat Commun 2020: Neural integration underlying 
%naturalistic prediction flexibly adapts to varying sensory input rate

% Scripts for correlation between neuromagnetic estimates of sensory 
% history integration (Kprime) and behavioral history dependence parameters
%(Fig. 6)

%% 0.1 Prepare input param_NASTD_MEGeters
%0.1.1 Specify paths
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/MEG/NASTD_MEG_Matlab/GitHub/')
paths_NASTD_MEG = NASTD_MEG_SetPaths;
addpath(genpath(paths_NASTD_MEG.ScriptsDir));

paths_NASTD_MEG.Current_analysis = paths_NASTD_MEG.Analysis.HisTrack.Correlation_KprimeBehav;
paths_NASTD_MEG.Current_outputdata =  ([paths_NASTD_MEG.Current_analysis 'Data/']);
paths_NASTD_MEG.Current_outputfig = [paths_NASTD_MEG.Current_analysis 'Figs/'];

%0.1.2 Specify subject param_NASTD_MEGeters
sub = 'S1'; %placeholder
NASTD_MEG_SubInfo %load subject info file (var: si)
param_NASTD_MEG.subs = si.sub_list;
clear sub si

%0.1.3 Specify analysis and plotting param_NASTD_MEGeters
param_NASTD_MEG.tonedur_text = {'0.15' '0.3' '0.6'};

param_NASTD_MEG.plot.plot = 1; %Plot Figures?
param_NASTD_MEG.plot.save = 0; %Save Figures?
param_NASTD_MEG.plot.pval = 0.025; %pval determining sign. sensors in topoplot (two-tailed test at p < 0.05)

param_NASTD_MEG.cluster.NumReps = 1000; %Number permutations for null distribution for cluster correction
param_NASTD_MEG.cluster.pval_clusterdef = 0.05; %pvalue used to define sensor-clusters

param_NASTD_MEG.KprimeComparison.ClusterSOI = {'PredictiveProcessingCluster'}; 

%% 1) Load in data
path_input_Kprime = [paths_NASTD_MEG.Analysis.HisTrack.Kprime_Computation 'SummaryStruct/'] ; 
path_input_Behav  = [paths_NASTD_MEG.Analysis.Behavior 'FTPLrating/Statistics/']; 

%Neural history tracking effect - Exp Kprime data 
%(SSub & GAvg Exp k-prime vals from summary struct)
load([path_input_Kprime 'KprimeSummaryStruct_SSubsGAvg.mat'])
clear Kprime_GAvg
%Average kprime vals across ToneDur for first 3 time windows
for i_win = 1:length(Kprime_AllSub.Exp{1})
    for i_sub = 1:size(Kprime_AllSub.Exp{1}{i_win},1)
        ExpKprime_AllSub_AvgToneDur{i_win}(i_sub,:) = ...
            mean([Kprime_AllSub.Exp{1}{i_win}(i_sub,:); ...
            Kprime_AllSub.Exp{2}{i_win}(i_sub,:); ...
            Kprime_AllSub.Exp{3}{i_win}(i_sub,:)]);
    end
end
    
%Behavioral history dependence effect(computed as F-statistic from 
%p34 x p*34 interaction effect from 3-way RM ANOVA; DV: FTPLratings)
load([path_input_Behav 'FTPLrating_Stats_Subs.mat']);
for i_subs = 1:length(param_NASTD_MEG.subs)
    Fstatinteraction_FTPLrating(i_subs,1) = ...
        FTPLstat.SingleSub.ANOVA3.FP_Inter_p34predp34{i_subs,1};
end
clear FTPLstat
 
%% 2) Correlate Kprime and behavioral parameter across subjects for each sensor
for i_win = 1:length(ExpKprime_AllSub_AvgToneDur)
    [rho_CorrexpKFstat_allSens{i_win}, pval_CorrexpKFstat_allSens{i_win}] = ...
        corr(ExpKprime_AllSub_AvgToneDur{i_win}, ...
        Fstatinteraction_FTPLrating(:,1), ...
        'type', 'Spearman','tail','both');
end

%% 3) Correlate Kprime and FTPLratings-Fvalues across subjects 
%for predetermined predictive processing clusters

%Load in file containing all sensor-clusters
load([paths_NASTD_MEG.Analysis.HisTrack.Kprime_Computation 'clusterSOI.mat'])

for i_tonedur = 1:length(param_NASTD_MEG.tonedur_text)

    filter_Subs = 1:length(param_NASTD_MEG.subs); %All subs

    for i_clustertype = 1:length(param_NASTD_MEG.KprimeComparison.ClusterSOI)
        clusterlabel = [param_NASTD_MEG.KprimeComparison.ClusterSOI{i_clustertype}];        
        for i_clusterwindow = 1:length(ClusterSOI.(clusterlabel))
            for i_clusternum = 1:length(ClusterSOI.(clusterlabel){i_clusterwindow})
                if ~isempty(ClusterSOI.(clusterlabel){i_clusterwindow}{i_clusternum})
                    
                    [rho_CorrexpKFstat_perSOI{i_clusterwindow,i_clusternum}, ...
                        pval_CorrexpKFstat_perSOI{i_clusterwindow,i_clusternum}] = ...
                        corr(mean(ExpKprime_AllSub_AvgToneDur{i_clusterwindow}...
                        (:,ClusterSOI.(clusterlabel){i_clusterwindow}{i_clusternum}'),2), ...
                        Fstatinteraction_FTPLrating(:,1), ...
                        'type', 'Spearman','tail','both');
                    
                end
            end
        end
    end
end

%% 4) Correct correlation p-values for multiple comparisons over time windows via FDR
%Combine p-values per TD (across TW and cluster)
allp_uncorr = [];
max_win = length(ExpKprime_AllSub_AvgToneDur); %check max. num of windows based on Kprime data
if max_win > length(ClusterSOI.(clusterlabel)) 
    %If number of windows based on Kprime data is larger than cluster-defined windows, restrict to latter
    max_win =  length(ClusterSOI.(clusterlabel));
end

counter_tonedur = 1;
counter_entries = 1;
for i_win = 1:max_win
    for i_cluster = 1:size(pval_CorrexpKFstat_perSOI,2)
        if ~isempty(pval_CorrexpKFstat_perSOI{i_win,i_cluster})
            allp_uncorr(counter_tonedur,1) = pval_CorrexpKFstat_perSOI{i_win,i_cluster};
            counter_tonedur = counter_tonedur +1;
        end
    end
end
allp_FDRcorr = mafdr(allp_uncorr,'BHFDR', true); %Note: function mafdr requires Matlab2017a or the Bioinformatics Toolbox
for i_win = 1:max_win
    for i_cluster = 1:size(pval_CorrexpKFstat_perSOI,2)
        if ~isempty(pval_CorrexpKFstat_perSOI{i_win,i_cluster})
            pvalFDRcorrected_CorrexpKFstat_perSOI{i_win,i_cluster} = ...
                allp_FDRcorr(counter_entries);
            counter_entries = counter_entries +1;
        end
    end
end

%% 5) Plot correlations for each predetermined sensor-cluster
%(Corresponds to Fig. 6a)
NASTD_MEG_Corr_PlotCorrSOI...
    (ExpKprime_AllSub_AvgToneDur, Fstatinteraction_FTPLrating, ... 
    ClusterSOI, ...
    rho_CorrexpKFstat_perSOI, pval_CorrexpKFstat_perSOI, pvalFDRcorrected_CorrexpKFstat_perSOI, ...
    param_NASTD_MEG, ...
    paths_NASTD_MEG);

%% 6) Compute cluster-corrected correlations for entire sensor array
%(Corresponds to Fig. 6b)
NASTD_MEG_Corr_ClusterCorr_Group...
    (Kprime_AllSub, Fstatinteraction_FTPLrating, ... 
    rho_CorrexpKFstat_allSens, pval_CorrexpKFstat_allSens, ...
    param_NASTD_MEG, ...
    paths_NASTD_MEG);