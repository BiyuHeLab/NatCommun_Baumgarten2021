%% Baumgarten et al. Nat Commun 2020: Neural integration underlying
%naturalistic prediction flexibly adapts to varying sensory input rate

% Scripts for analysis of ERF effects (neuromagnetic activity over the
% course of the tone sequence presentation as a functionof p*34 for
% predictive processing clusters and early sensory filters)
%(Fig. 3b, Supplementary Figure 3)

%% 0.1 Prepare input parameters
%0.1.1 Specify paths
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/MEG/NASTD_MEG_Matlab/GitHub/') %Set base dir
paths_NASTD_MEG = NASTD_MEG_SetPaths;
addpath(genpath(paths_NASTD_MEG.ScriptsDir));

paths_NASTD_MEG.Current_inputdata = paths_NASTD_MEG.Rawdata.MEG;
paths_NASTD_MEG.Current_analysis = paths_NASTD_MEG.Analysis.ERF;
paths_NASTD_MEG.Current_outputdata =  ([paths_NASTD_MEG.Current_analysis]);
paths_NASTD_MEG.Current_outputfig = [paths_NASTD_MEG.Current_outputdata 'Figs/'];

%0.1.2 Specify subject parameters
sub = 'S1'; %placeholder
NASTD_MEG_SubInfo %load subject info file (var: si)
param_NASTD_MEG.subs = si.sub_list;
clear sub si

%0.1.3 Specify analysis and plotting parameters
param_NASTD_MEG.tonedur_text = {'0.15' '0.3' '0.6'};
param_NASTD_MEG.samplefreq = 600;
param_NASTD_MEG.plot.plot = 1; %Plot Figures?
param_NASTD_MEG.plot.save = 0; %Save Figures?

param_NASTD_MEG.M100.TW = [0.075 0.125]; %in S
param_NASTD_MEG.M100.IncludedTones = [1:34]; %All tones
param_NASTD_MEG.M100.Subs{1} = [1 3 4 5 6 7 8 9 10 11 13 14 15 16 17 18 19 20]; %Subs that show reliable M100 response for 150ms TD
param_NASTD_MEG.M100.Subs{2} = [1 2 3 4 5 6 7 8 9 10 12 13 14 15 16 17 18 19 20];%Subs that show reliable M100 response for 300ms TD
param_NASTD_MEG.M100.Subs{3} = [1 2 3 5 6 7 8 9 10 12 13 14 15 16 17 18 19 20];%Subs that show reliable M100 response for 600ms TD

param_NASTD_MEG.stat.NumReps = 1000; %Number permutations for null distribution for cluster correction
param_NASTD_MEG.stat.pval_clusterstat = 0.05; %pvalue used for cluster statistics

%Sensor Cluster used for sensors of interest (SOI) selection
param_NASTD_MEG.SOI.Type = {'PredictiveProcessingCluster'};
param_NASTD_MEG.SOI.pval_load = 0.05; %pval used to determine which cluster-files are loaded
param_NASTD_MEG.SOI.nPerm = 1000; %Number of permutations used for cluster computation - determines cluster data load in
param_NASTD_MEG.SOI.pval_select = 0.025; %pval used to determine which clusters are selected from loaded file
param_NASTD_MEG.SOI.windows = 1:3;

%% 1. Compute and plot ERFs for Sensors of Interest (SOI, Predictive Processing Cluster)
%1.1 Determine Input: Sensors of Interest (SOI) PredictiveProcessingCluster per Tone Duration
for i_tonedur = 1:length(param_NASTD_MEG.tonedur_text)
    
    clusterSOI = ...
        NASTD_MEG_ERF_GetSensorIndexfromCluster...
        (param_NASTD_MEG.tonedur_text{i_tonedur}, param_NASTD_MEG, paths_NASTD_MEG);
    
    for i_win_cluster = 1:length(clusterSOI)
        for i_cluster = 1:length(clusterSOI{i_win_cluster}.i_SensperCluster)
            index_clusterSOI{i_tonedur}{i_win_cluster}{i_cluster} = ...
                clusterSOI{i_win_cluster}.i_SensperCluster{i_cluster};
        end
    end
end

for i_tonedur = 1:length(param_NASTD_MEG.tonedur_text)
    %1.2) Distinguish and read out trials per Tone Duration & subject based ...
    %on p*34
   NASTD_MEG_ERF_SelectTrials...
        (i_tonedur, param_NASTD_MEG, paths_NASTD_MEG);
    
    %1.3 Compute ERF for p*34 trial and sensors-subselection for each subject and place
    %them in common-FT-style struct
    NASTD_MEG_ERF_ComputeERFperCondition...
        (i_tonedur, ...
        index_clusterSOI, ...
        param_NASTD_MEG, paths_NASTD_MEG);
    
    %1.4 Perform group-level statistical comparison for each time point (FT-cluster corrected)
    %Between low vs. high p*34 trials, focused on and averaged across SOI
    ERF_GroupStat_LowvsHighPredp34{i_tonedur} = ...
        NASTD_MEG_ERF_ERFGroupStat...
        (i_tonedur, ...
        param_NASTD_MEG, paths_NASTD_MEG);
    
    %1.5 Plot ERFs over the course of the tone sequence and highlight
    %samples showing sign. differences between low vs. high p*34
    NASTD_MEG_ERF_PlotERF_GroupStat...
        (ERF_GroupStat_LowvsHighPredp34{i_tonedur}, ...
        i_tonedur, ...
        param_NASTD_MEG,paths_NASTD_MEG)
end

%% 2. Compute and plot ERFs for M100 Spatial Filter
for i_tonedur = 1:length(param_NASTD_MEG.tonedur_text)
    
    %2.1 Compute M100 spatial filter weights for each subject and 
    %plot summary figure to determine subjects with valid M100 response
    NASTD_MEG_ERF_CompM100SpatialFilter...
        (i_tonedur, param_NASTD_MEG, paths_NASTD_MEG);

    %2.2 Perform group-level statistical comparison for weighted 
    %low vs. high p*34 trials 
    NASTD_MEG_ERF_StatM100SpatialFilter...
        (i_tonedur, ...
        param_NASTD_MEG, paths_NASTD_MEG);
end
