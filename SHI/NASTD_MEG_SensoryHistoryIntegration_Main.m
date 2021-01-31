%% Baumgarten et al. Nat Commun 2020: Neural integration underlying 
%naturalistic prediction flexibly adapts to varying sensory input rate

% Scripts for analysis of neural sensory history (SHI) effects 
%(corresponds to Fig. 5, Supplementary Figure 4, Supplementary Figure 5, 
%Supplementary Figure 6, Supplementary Figure 7c)

%% 0.1 Prepare input parameters
%0.1.1 Specify paths
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/MEG/NASTD_MEG_Matlab/GitHub/')
paths_NASTD_MEG = NASTD_MEG_SetPaths;
addpath(genpath(paths_NASTD_MEG.ScriptsDir));

paths_NASTD_MEG.Current_inputdata = paths_NASTD_MEG.Rawdata.MEG;
paths_NASTD_MEG.Current_analysis = paths_NASTD_MEG.Analysis.HisTrack.Kprime_Computation;
paths_NASTD_MEG.Current_outputdata =  ([paths_NASTD_MEG.Current_analysis]);
paths_NASTD_MEG.Current_outputfig = [paths_NASTD_MEG.Current_outputdata 'Figs/'];

%0.1.2 Specify subject param_NASTD_MEGeters
sub = 'S1'; %placeholder
NASTD_MEG_SubInfo %load subject info file (var: si)
param_NASTD_MEG.subs = si.sub_list;
clear sub si

%0.1.3 Specify analysis and plotting param_NASTD_MEGeters
param_NASTD_MEG.tonedur_text = {'0.15' '0.3' '0.6'};
param_NASTD_MEG.betalevel_text = {'0.5','0.99','1.5'};

param_NASTD_MEG.plot.plot = 1; %Plot Figures?
param_NASTD_MEG.plot.save = 0; %Save Figures?
param_NASTD_MEG.plot.pval = 0.05; %pval determining sign. sensors in topoplot (one-tailed test at p < 0.05)

param_NASTD_MEG.NumFolds = 6; %Number of folds used for cross-validation
param_NASTD_MEG.NumReps = 100; %Number of repetitions to create shuffled null distribution

param_NASTD_MEG.cluster.NumReps = 1000; %Number permutations used for cluster correction
param_NASTD_MEG.cluster.pval_clusterdef = 0.05; %p-value used to define sign. sensor-clusters

param_NASTD_MEG.KprimeComparison.coord_StartPoint = [0,0,0]; %Start point
param_NASTD_MEG.KprimeComparison.coord_EndPoint_DurationLine = [10 5 2.5];%End point Duration line
param_NASTD_MEG.KprimeComparison.coord_EndPoint_InformationLine = [10 10 10];%End point Information line
param_NASTD_MEG.KprimeComparison.NumSamples4ShuffGroupDist = 1000; 
param_NASTD_MEG.KprimeComparison.ClusterSOI = ...
    {'PredictiveProcessingCluster',...
    'SensorCluster_DurationLineEffect',...
    'SensorCluster_InformationLineEffect'}; 
param_NASTD_MEG.KprimeComparison.AnalysisParam = ...
    {'Norm','AngleInformationLine','AngleDurationLine'}; 

%% 1) Compute, save, and optionally plot Kprime effects for experimental data
%1.1) Compute k-values (number of past tones that explains current 
%neuromagnetic activity) for experimental data for each tone duration
for i_sub = 1:length(param_NASTD_MEG.subs)
    for i_tonedur = 1:length(param_NASTD_MEG.tonedur_text)
        NASTD_MEG_SHI_CompExpK...
            (param_NASTD_MEG.subs{i_sub}, ...
            param_NASTD_MEG.tonedur_text{i_tonedur},...
            param_NASTD_MEG.NumFolds, ...
            paths_NASTD_MEG);    
    end
end

%1.2 Combine k-values across folds to get k-prime
for i_sub = 1:length(param_NASTD_MEG.subs)
    for i_tonedur = 1:length(param_NASTD_MEG.tonedur_text)
        NASTD_MEG_SHI_CompExpKprime...
            (param_NASTD_MEG.subs{i_sub}, ...
            param_NASTD_MEG.tonedur_text{i_tonedur},...
            param_NASTD_MEG.plot.plot, param_NASTD_MEG.plot.save,...
            paths_NASTD_MEG);
    end
end
        
%1.3 Average k-prime values across subjects and plot resulting topo
for i_tonedur = 1:length(param_NASTD_MEG.tonedur_text)
    NASTD_MEG_SHI_GroupEXPKprime...
        (param_NASTD_MEG.subs, ...
        param_NASTD_MEG.tonedur_text{i_tonedur}, ...
        param_NASTD_MEG.plot.plot, param_NASTD_MEG.plot.save,...
        paths_NASTD_MEG);
end

%% 2) Compute and save Kprime null distributions
%2.1 Compute k-values (number of past tones that explains current 
%neuromagnetic activity) for shuffled data for each tone duration
for i_sub = 1:length(param_NASTD_MEG.subs)
    NASTD_MEG_SHI_CompShuffK...
        (param_NASTD_MEG.subs{i_sub}, ...
        param_NASTD_MEG.tonedur_text,...
        param_NASTD_MEG.NumFolds, param_NASTD_MEG.NumReps,...
        paths_NASTD_MEG);
end

%% 3) Determine sign. Kprime values by comparing experimental vs. null Kprime values
%3.1 Combine shuffled folds to get shuffled k-prime values. Compare Exp +
%Shuffled k-prime to get single-subject stat. assessment of k-prime values
for i_sub = 1:length(param_NASTD_MEG.subs)
    for i_tonedur = 1:length(param_NASTD_MEG.tonedur_text)
        NASTD_MEG_SHI_CompareExpShuffKprime...
            (param_NASTD_MEG.subs{i_sub}, ...
            param_NASTD_MEG.tonedur_text{i_tonedur},...
            param_NASTD_MEG.plot.plot, param_NASTD_MEG.plot.pval, param_NASTD_MEG.plot.save,...
            paths_NASTD_MEG);       
    end
end

%3.2 Create Group-level k-prime null distribution, compute and plot 
%Cluster-corrected statistics as topoplot (corresponds to Supplementary Fig. 4)
for i_tonedur = 1:length(param_NASTD_MEG.tonedur_text)
    NASTD_MEG_SHI_ClusterCorrKprime_Group...
        (param_NASTD_MEG.subs, param_NASTD_MEG.tonedur_text{i_tonedur},...
        param_NASTD_MEG.cluster.NumReps, param_NASTD_MEG.cluster.pval_clusterdef, ...
        param_NASTD_MEG.plot.plot, param_NASTD_MEG.plot.pval, param_NASTD_MEG.plot.save,...
        paths_NASTD_MEG);
end

%3.3 Create summary structure containing all single-subject Kprime
%parameters (exp & shuffled) for all TD
NASTD_MEG_SHI_KprimeSummaryStruct...
    (param_NASTD_MEG.subs, param_NASTD_MEG.tonedur_text, ...
    paths_NASTD_MEG);

%% 4) Compare Kprime values across tone duration conditions
%4.1 Combine Kprime across all tone durations and compute
%vector norm and angle, and significance per sensor and sensor cluster
NASTD_MEG_SHI_CompKprimeNormAngle_Group...
    (param_NASTD_MEG, ...
    paths_NASTD_MEG);

%4.2 Plot ouput of statistical exp vs. shuff comparison for exp vs. shuff 
%vector norm and angle (corresponds to Fig. 5 and Supplementary Fig. 5)
NASTD_MEG_SHI_PlotKprimeNormAngle_Group...
    (param_NASTD_MEG, ...
    paths_NASTD_MEG);
    
%4.3 Perform cluster-correction for exp vs. shuff data across the entire 
%sensor array (corresponds to Fig. 5 and Supplementary Fig. 6)
NASTD_MEG_SHI_ClusterCorrKprimeNormAngle_Group...
    (param_NASTD_MEG, ...
    paths_NASTD_MEG);

%4.4 Project Kprime values in 2D space and plot them in 2D coordinate system
%(corresponds to Fig. 5 and Supplementary Fig. 5)
NASTD_MEG_SHI_CompKprimeNormAngle_2Dspace_Group...
    (param_NASTD_MEG, ...
    paths_NASTD_MEG);

%% 5) Compute Kprime effects based on sequence beta (not tone duration)
%5.1) Compute and combine experimental k-values (neural activity as 
%function of prior tone pitch) or experimental data for each sequence beta level
for i_sub = 1:length(subs) %Loop subjects
    for i_betalevel = 1:length(param_NASTD_MEG.betalevel_text) 
        NASTD_MEG_SHI_CompExpK_perBeta...
            (param_NASTD_MEG.subs{i_sub}, param_NASTD_MEG.betalevel_text{i_betalevel}, ...
            param_NASTD_MEG.tonedur_text,...
            param_NASTD_MEG.NumFolds, ...
            paths_NASTD_MEG); 

    end
end

%5.2) Compute shuffled k-values for each sequence beta level
for i_sub = 1:length(subs) %Loop subjects
    for i_betalevel = 1:length(param_NASTD_MEG.betalevel_text) 
        NASTD_MEG_SHI_CompShuffK_perBeta...
            (param_NASTD_MEG.subs{i_sub}, param_NASTD_MEG.betalevel_text{i_betalevel}, ...
            param_NASTD_MEG.tonedur_text, ...
            param_NASTD_MEG.NumFolds, param_NASTD_MEG.NumReps, ...
            paths_NASTD_MEG); 

    end
end

%5.3 Compare Exp + Shuffled k-prime values to get single-subject stat. 
%assessment of k-prime values per beta level 
for i_sub = 1:length(subs) %Loop subjects
    for i_betalevel = 1:length(param_NASTD_MEG.betalevel_text) 
        NASTD_MEG_SHI_CompareExpShuffKprime_perBeta...
            (param_NASTD_MEG.subs{i_sub}, param_NASTD_MEG.betalevel_text{i_betalevel}, ...
            param_NASTD_MEG.tonedur_text, ...
            param_NASTD_MEG.NumFolds, param_NASTD_MEG.NumReps, ...
            paths_NASTD_MEG); 

    end
end

%5.4 Create Group-level k-prime null distribution, Compute and plot 
%Cluster-corrected statistics as topoplot
%(corresponds to Supplementary Fig. 7)
for i_betalevel = 1:length(param_NASTD_MEG.betalevel_text) 
    NASTD_MEG_SHI_ClusterCorrKprime_Group_perBeta...
        (param_NASTD_MEG.subs, param_NASTD_MEG.betalevel_text{i_betalevel}, ...
        param_NASTD_MEG.tonedur_text,...
        param_NASTD_MEG.cluster.NumReps, param_NASTD_MEG.cluster.pval_clusterdef, ...
        param_NASTD_MEG.plot.pval, param_NASTD_MEG.plot.save,...
        paths_NASTD_MEG);
end