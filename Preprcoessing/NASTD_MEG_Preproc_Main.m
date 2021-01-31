%% Baumgarten et al. Nat Commun 2020: Neural integration underlying 
%naturalistic prediction flexibly adapts to varying sensory input rate

% Scripts for preprocessing of MEG data

%% 0.1 Prepare input parameters
%0.1.1 Specify paths
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/MEG/NASTD_MEG_Matlab/GitHub/') %Set base dir
paths_NASTD_MEG = NASTD_MEG_SetPaths;
addpath(genpath(paths_NASTD_MEG.ScriptsDir));

paths_NASTD_MEG.Current_inputdata   = paths_NASTD_MEG.Rawdata.MEG;
paths_NASTD_MEG.Current_outputdata  = paths_NASTD_MEG.Rawdata.MEG;
paths_NASTD_MEG.Current_outputfig   = [paths_NASTD_MEG.Current_outputdata 'Figs/'];

%0.1.2 Specify subject parameters
sub = 'S1'; %placeholder
NASTD_MEG_SubInfo %load subject info file (var: si)
param_NASTD_MEG.subs = si.sub_list;
clear sub si

%% 1) Read out events, adjust trigger delay to real stimulus onset, and save adjusted events
for i_sub = 1:length(param_NASTD_MEG.subs)
    
    NASTD_MEG_Preproc_AdjustStimOnset...
        (param_NASTD_MEG.subs{i_sub}, paths_NASTD_MEG)
    
    NASTD_MEG_Preproc_SaveEvents...
        (param_NASTD_MEG.subs{i_sub}, paths_NASTD_MEG)   
    
%% 2) Visualize raw data and check for bad sensors

    NASTD_MEG_Preproc_CheckRawData...
        (param_NASTD_MEG.subs{i_sub}, paths_NASTD_MEG)

      
%% 3) Perform ICA to remove EOG/cardiac/technical artifacts

    NASTD_MEG_Preproc_ICAcalculation...
        (param_NASTD_MEG.subs{i_sub}, paths_NASTD_MEG)
    
    NASTD_MEG_Preproc_ICAselection...
        (param_NASTD_MEG.subs{i_sub}, paths_NASTD_MEG)
    
    NASTD_MEG_Preproc_ICAremove...
        (param_NASTD_MEG.subs{i_sub}, paths_NASTD_MEG)  
    
end