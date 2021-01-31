%% Baumgarten et al. Nat Commun 2020: Neural integration underlying 
%naturalistic prediction flexibly adapts to varying sensory input rate

% Scripts for analysis of neural prediction effects 
%(Fig. 3a, Supplementary Figure 1, 2, 7)

%% 0.1 Prepare input parameters
%0.1.1 Specify paths
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/MEG/NASTD_MEG_Matlab/GitHub/') %Set base dir
paths_NASTD_MEG = NASTD_MEG_SetPaths;
addpath(genpath(paths_NASTD_MEG.ScriptsDir));

paths_NASTD_MEG.Current_inputdata = paths_NASTD_MEG.Rawdata.MEG;
paths_NASTD_MEG.Current_analysis = paths_NASTD_MEG.Analysis.Prediction;
paths_NASTD_MEG.Current_outputdata =  ([paths_NASTD_MEG.Current_analysis]);
paths_NASTD_MEG.Current_outputfig = [paths_NASTD_MEG.Current_outputdata 'Figs/'];

%0.1.2 Specify subject parameters
sub = 'S1'; %placeholder
NASTD_MEG_SubInfo %load subject info file (var: si)
param_NASTD_MEG.subs = si.sub_list;
clear sub si

%0.1.3 Specify analysis and plotting parameters
param_NASTD_MEG.tonedur_text = {'0.15' '0.3' '0.6'};
param_NASTD_MEG.betalevel_text = {'0.5','0.99','1.5'};

param_NASTD_MEG.plot.plot = 1; %Plot Figures?
param_NASTD_MEG.plot.save = 0; %Save Figures?
param_NASTD_MEG.plot.pval = 0.025; %pval determining sign. sensors in topoplot (two-tailed test at p < 0.05)

param_NASTD_MEG.predictive_sequencerange = [33]; %final tone defining range for tone pitches based on which prediction should be computed (i.e., 1-33)
param_NASTD_MEG.toneIndex = 33; %tone index indicating for which tone ERF/ERP should be computed, which is then used in regression analysis of prediciton effect

param_NASTD_MEG.cluster.NumReps = 1000; %Number permutations for null distribution for cluster correction
param_NASTD_MEG.cluster.pval_clusterdef = 0.05; %pvalue used to define sensor-clusters

%% 1) Compute and save prediction effects (optionally plot results)
%1.1) Compute prediction effect (neural activity at tone 33 that is
%predictive of math. expected tone/p*34) for each subject and tone duration
for i_sub = 1:length(param_NASTD_MEG.subs)
    for i_tonedur = 1:length(param_NASTD_MEG.tonedur_text) 
        NASTD_MEG_Pred_CompPred_Subs...
            (param_NASTD_MEG.subs{i_sub}, param_NASTD_MEG.tonedur_text{i_tonedur},...
            param_NASTD_MEG.predictive_sequencerange, param_NASTD_MEG.toneIndex, ...
            param_NASTD_MEG.plot.plot, param_NASTD_MEG.plot.pval, param_NASTD_MEG.plot.save,...
            paths_NASTD_MEG);
    end
end

%1.1B) Alternatively, compute prediction effect for restricted trial selection
param_NASTD_MEG.selected_trials = 'firsthalf'; %firsthalf, lasthalf
for i_sub = 1:length(param_NASTD_MEG.subs)
    for i_tonedur = 1:length(param_NASTD_MEG.tonedur_text)
        NASTD_MEG_Pred_CompPred_TrialSelection_Subs...
            (param_NASTD_MEG.subs{i_sub}, param_NASTD_MEG.tonedur_text{i_tonedur},...
            param_NASTD_MEG.predictive_sequencerange, param_NASTD_MEG.toneIndex, ...
            param_NASTD_MEG.selected_trials,...
            param_NASTD_MEG.plot.plot, param_NASTD_MEG.plot.pval, param_NASTD_MEG.plot.save,...
            paths_NASTD_MEG);
    end
end

%1.1C) Alternatively, compute prediction effect per sequence beta level
for i_sub = 1:length(param_NASTD_MEG.subs)
    for i_betalevel = 1:length(param_NASTD_MEG.betalevel_text)         
        NASTD_MEG_Pred_CompPred_Beta_Subs ...
            (param_NASTD_MEG.subs{i_sub}, param_NASTD_MEG.betalevel_text{i_betalevel}, ...
            param_NASTD_MEG.tonedur_text, ...
            param_NASTD_MEG.predictive_sequencerange, param_NASTD_MEG.toneIndex, ...
            paths_NASTD_MEG);        
    end
end

%1.2) Compute cluster-corrected single-subject effect (to get individual prediction sensor-clusters)
for i_sub = 1:length(param_NASTD_MEG.subs)
    for i_tonedur = 1:length(param_NASTD_MEG.tonedur_text)
       NASTD_MEG_Pred_ClusterCorrPred_Subs...
            (param_NASTD_MEG.subs{i_sub},  param_NASTD_MEG.tonedur_text{i_tonedur},...
            param_NASTD_MEG.cluster.pval_clusterdef, param_NASTD_MEG.cluster.NumReps, ...
            param_NASTD_MEG.plot.plot, param_NASTD_MEG.plot.pval, param_NASTD_MEG.plot.save,...
            paths_NASTD_MEG);
    end
end

%1.3) Compute and plot cluster-corrected group-level effect 
%(Corresponds to Fig. 3a middle inset, Supplementary Fig. 1, 2)
for i_tonedur = 1:length(param_NASTD_MEG.tonedur_text)   
    NASTD_MEG_Pred_ClusterCorrPred_Group...
            (param_NASTD_MEG.subs, param_NASTD_MEG.tonedur_text{i_tonedur},...
            param_NASTD_MEG.cluster.pval_clusterdef, param_NASTD_MEG.cluster.NumReps, ...
            param_NASTD_MEG.plot.plot, param_NASTD_MEG.plot.pval, param_NASTD_MEG.plot.save,...
            paths_NASTD_MEG);    
end

%1.3C) Compute and plot cluster-corrected group-level effect for trials
%based on sequence beta level subselection 
%(Corresponds to Supplementary Fig. 7b)
for i_betalevel = 1:length(param_NASTD_MEG.betalevel_text)
    NASTD_MEG_Pred_ClusterCorrPred_Beta_Group...
            (param_NASTD_MEG.subs,  param_NASTD_MEG.betalevel_text{i_betalevel}, ...
            param_NASTD_MEG.tonedur_text, ...
            param_NASTD_MEG.cluster.pval_clusterdef, param_NASTD_MEG.cluster.NumReps, ...
            param_NASTD_MEG.plot.plot, param_NASTD_MEG.plot.pval, param_NASTD_MEG.plot.save,...
            paths_NASTD_MEG);     
end

%% 2) Compute timelocked activity during p33 as a function of p*34
%2.1 Compute and plot timelocked activity during p33 as a function of p*34
%Corresponds to Fig. 3a right inset
for i_tonedur = 1:length(param_NASTD_MEG.tonedur_text)   
    NASTD_MEG_Pred_ERFp33perpredp34_Group...
            (param_NASTD_MEG.subs,  param_NASTD_MEG.tonedur_text{i_tonedur},...
            param_NASTD_MEG.toneIndex, ...
            param_NASTD_MEG.plot.plot, param_NASTD_MEG.plot.save,...
            paths_NASTD_MEG);    
end