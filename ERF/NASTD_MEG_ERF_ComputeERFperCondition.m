function NASTD_MEG_ERF_ComputeERFperCondition...
    (i_tonedur, index_clusterSOI, param, paths_NASTD_MEG)

%Aim: Compute timelocked event-related fields (ERF) for neural data from a
%specific sensor-selection (predictive processing clusters) and a specific
%trial-subselection (p*34 low vs. high). Combine ERFs from all subjects in
%a common output struct.
path_input = [paths_NASTD_MEG.Current_outputdata 'Data/Ssub/' param.tonedur_text{i_tonedur} 'sTD/'];
path_output = [paths_NASTD_MEG.Current_outputdata 'Data_PredProcCluster/Ssub/' param.tonedur_text{i_tonedur} 'sTD/'];
mkdir(path_output)

for i_sub = 1:length(param.subs)
    
    ERF_avgTrialsperCond = struct;

    disp(['Sub' num2str(param.subs{i_sub}) ' - TD: ' param.tonedur_text{i_tonedur} ' - Start processing ERF'])
    
    %Load in single trial data
    load([path_input param.subs{i_sub} '_MEGsingletrials_' ...
        param.tonedur_text{i_tonedur} 'sTD.mat']);
    
    %Determine tone-start-samples
    nSamplesPerTone = str2double(param.tonedur_text{i_tonedur}) * param.samplefreq;
    index_T1start = find(MEGactivity_perTrial.time{1} == 0); %determine start of first tone
    Samples_AllToneStart = [index_T1start : nSamplesPerTone : (index_T1start + (34*nSamplesPerTone))];
    
    %Specify SOI and average activity across these sensors
    if strcmp(param.tonedur_text{i_tonedur},'0.15') %For 150ms TD
        index_SOI = index_clusterSOI{i_tonedur}{3}{1}; %Take TW3, right cluster
    elseif strcmp(param.tonedur_text{i_tonedur},'0.3') %For 300ms TD
        index_SOI = index_clusterSOI{i_tonedur}{2}{2}; %Take TW2, right cluster
    elseif strcmp(param.tonedur_text{i_tonedur},'0.6') %For 600ms TD
        index_SOI = index_clusterSOI{i_tonedur}{3}{1}; %Take TW3, right cluster
    end
    
    cfg                    = [];
    cfg.channel            = index_SOI;
    cfg.avgoverchan        = 'yes';
    MEGactivity_perTrial_AvgSOI{i_sub} = ...
        ft_selectdata(cfg, MEGactivity_perTrial);
    
    %Select p*34 (predp34) sub-conditions
    conditions = fieldnames(MEGactivity_perTrial_AvgSOI{i_sub}.trialindices);
    for i_cond = 1:length(conditions)
        cfg                             = [];
        cfg.trials                      = MEGactivity_perTrial_AvgSOI{i_sub}.trialindices.(conditions{i_cond});
        MEGactivity_perTrialperCond_AvgSOI{i_sub}.(conditions{i_cond}) = ...
            ft_preprocessing(cfg, MEGactivity_perTrial_AvgSOI{i_sub});
                      
        %Baseline Correct for prestimulus baseline (-0.5 to 0s0
        TOI_baseline_samples = [1 Samples_AllToneStart(1)];
        TOI_baseline_secs = ...
            [MEGactivity_perTrialperCond_AvgSOI{i_sub}.(conditions{i_cond}).time{1}(TOI_baseline_samples(1)) ...
            MEGactivity_perTrialperCond_AvgSOI{i_sub}.(conditions{i_cond}).time{1}(TOI_baseline_samples(2))];
                
        cfg                     = [];
        cfg.demean              = 'yes';
        cfg.baselinewindow      = TOI_baseline_secs;
        MEGactivity_perTrialperCond_AvgSOI{i_sub}.(conditions{i_cond}) = ...
            ft_preprocessing(cfg,MEGactivity_perTrialperCond_AvgSOI{i_sub}.(conditions{i_cond}));
        
        %2.8 Compute ERF for condition-selected trials
        cfg                     = [];
        cfg.keeptrials          = 'no';
        ERF_avgTrialsperCond.(conditions{i_cond}) = ...
            ft_timelockanalysis(cfg, MEGactivity_perTrialperCond_AvgSOI{i_sub}.(conditions{i_cond}));
        
    end
    
    %Save single-trial ERFs
    savefile = [path_output  param.subs{i_sub} '_PredProcCluster_' ...
        param.tonedur_text{i_tonedur} 'sTD.mat'];
    save(savefile, 'ERF_avgTrialsperCond', '-v7.3');  
    
    clear ERF_avgTrialsperCond MEGactivity_perTrial
end

disp(['TD: ' param.tonedur_text{i_tonedur} ' - Finished processing All subs'])
fprintf('\n\n');

end