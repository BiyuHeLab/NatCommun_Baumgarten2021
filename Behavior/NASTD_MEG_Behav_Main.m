%% Baumgarten et al. Nat Commun 2020: Neural integration underlying 
%naturalistic prediction flexibly adapts to varying sensory input rate

% Scripts for analysis of behavioral responses 
%(Fig. 2, Supplementary Figure 7a, Input for Fig. 6)

%% 0 Set input parameters
%Set paths
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/MEG/NASTD_MEG_Matlab/GitHub/') %Set base dir
paths_NASTD_MEG = NASTD_MEG_SetPaths;
addpath(genpath(paths_NASTD_MEG.ScriptsDir));

paths_NASTD_MEG.Current_inputdata = paths_NASTD_MEG.Rawdata.Behav;
paths_NASTD_MEG.Current_analysis = paths_NASTD_MEG.Analysis.Behavior;
paths_NASTD_MEG.Current_outputdata =  ([paths_NASTD_MEG.Current_analysis 'FTPLrating/']);
paths_NASTD_MEG.Current_outputfig = [paths_NASTD_MEG.Current_outputdata 'Figs/'];

%Specify subject parameters
sub = 'S1'; %placeholder
NASTD_MEG_SubInfo %load subject info file (var: si)
param_NASTD_MEG.subs = si.sub_list;
clear sub si

%0.1.3 Specify analysis and plotting parameters
param_NASTD_MEG.tonedur_text = {'all' '0.15' '0.3' '0.6'};
param_NASTD_MEG.plot.plot = 1; %Plot Figures?
param_NASTD_MEG.plot.save = 0; %Save Figures?

%% 1. Plot and analyze Final Tone Pitch Likelihood (FTPL) ratings
%1.1 Single-subject level
for i_sub = 1:length(param_NASTD_MEG.subs)
    for i_tonedur = 1:length(param_NASTD_MEG.tonedur_text)
        NASTD_MEG_Behav_FTPL_Subs...
            (param_NASTD_MEG.subs{i_sub}, param_NASTD_MEG.tonedur_text{i_tonedur}, ...
            param_NASTD_MEG.plot.plot, param_NASTD_MEG.plot.save, ...
            paths_NASTD_MEG); 
    end
end

%1.2 Group-level 
%(corresponds to Fig. 2)
for i_tonedur = 1:length(param_NASTD_MEG.tonedur_text)
    NASTD_MEG_Behav_FTPL_Group...
        (param_NASTD_MEG.subs, param_NASTD_MEG.tonedur_text{i_tonedur}, ...
        param_NASTD_MEG.plot.plot, param_NASTD_MEG.plot.save, ...
            paths_NASTD_MEG); 
end

%% 2. Compute FTPL statistics
%2.1 Single-subject level
NASTD_MEG_Behav_FTPLStat_Subs...
    (param_NASTD_MEG.subs, param_NASTD_MEG.tonedur_text, ...
    paths_NASTD_MEG);

%2.2 Group-level
%Includes option to restrict trial-selection to first or last half of experiment
NASTD_MEG_Behav_FTPLStat_Group...
    (param_NASTD_MEG.subs, param_NASTD_MEG.tonedur_text, ...
    paths_NASTD_MEG);

%% 3. Compute FTPL statistics and plot group-level FTPL ratings per sequence beta 
%(corresponds to Supplementary Fig. 7a)
NASTD_MEG_Behav_FTPLperSeqBeta_Group...
    (param_NASTD_MEG.subs, param_NASTD_MEG.tonedur_text, ...
    paths_NASTD_MEG);
