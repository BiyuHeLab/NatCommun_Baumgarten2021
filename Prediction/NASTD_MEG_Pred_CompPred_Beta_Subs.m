function NASTD_MEG_Pred_CompPred_Beta_Subs...
    (sub, betalevel_input, tonedur_text, ...
    predictive_sequencerange, toneIndex, ...
    paths_NASTD_MEG)

%Aim: Compute and save neural prediction effect (How is neural activity at tone 33
%modulated by the expected value of tone 34, which itself depends on the
%previous tone sequence input?) for a specific beta level (across tone durations)

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath(genpath(paths_NASTD_MEG.ScriptsDir));
NASTD_MEG_SubInfo

%Data input
loadfile_MEG = [si.path_ica  'data_clean.mat'];%path to single-subject MEG data
loadfile_behav = si.path_behav;%path to single-subject behavioral/stimulus data

%Data output
path_save = [paths_NASTD_MEG.Current_outputdata 'PredEffect/SeqBeta' betalevel_input '/' sub '/'];

%% 1) Load and prepare input data
%Load preprocessed MEG data
disp('loading...')
load(loadfile_MEG); 
load(loadfile_behav); 

store_data_AllTD = data; %save copy of behav data

%Extract clean trials per tonedur, than append different tonedur
for i_tonedur = 1:length(tonedur_text)
    data = store_data_AllTD;
    data_MEG_perTD = extract_trials(data_clean, sub, str2double(tonedur_text{i_tonedur}));
    
    %Filter behavioral data to exclude trials that don't have MEG recording
    numoftrials = num2str(length(data.trialNum));
    
    dat_fields  = fieldnames(data);
    for i = 1:length(dat_fields)
        if eval(['length(data.' dat_fields{i} ') == 324'])
            eval(['data.' dat_fields{i} ' = data.' dat_fields{i} '(si.good_trials);']);
        end
    end
    
    stim_fields = fieldnames(data.stim);
    for i = 1:length(stim_fields)
        if eval(['length(data.stim.' stim_fields{i} ') == 324'])
            eval(['data.stim.' stim_fields{i} ' = data.stim.' stim_fields{i} '(si.good_trials);']);
        end
    end
    
    %Filter data_MEG for subjects where blocks have to be taken out
    dat_MEG_fields  = fieldnames(data_MEG_perTD);
    for i = 1:length(dat_MEG_fields)
        if eval(['length(data_MEG_perTD.' dat_MEG_fields{i} ') ==  ' numoftrials]) %360
            eval(['data_MEG_perTD.' dat_MEG_fields{i} ' = data_MEG_perTD.' dat_MEG_fields{i} '(si.good_trials);']);
        end
    end
    
    %Filter trials with selected ToneDur
    numoftrials = num2str(length(data.trialNum));
    
    f = data.stim.toneDur == str2double(tonedur_text{i_tonedur}); %filter to select trials for certain tone dur
    
    dat_fields  = fieldnames(data);
    for i = 1:length(dat_fields)
        if eval(['length(data.' dat_fields{i} ') ==' numoftrials])
            eval(['data.' dat_fields{i} ' = data.' dat_fields{i} '(f);']);
        end
    end
    
    stim_fields = fieldnames(data.stim);
    for i = 1:length(stim_fields)
        if eval(['length(data.stim.' stim_fields{i} ') ==' numoftrials])
            eval(['data.stim.' stim_fields{i} ' = data.stim.' stim_fields{i} '(f);']);
        end
    end
        
    MEG_fields = fieldnames(data_MEG_perTD);
    for i = 1:length(MEG_fields)
        if eval(['length(data_MEG_perTD.' MEG_fields{i} ') ==' numoftrials])
            eval(['data_MEG_perTD.' MEG_fields{i} ' = data_MEG_perTD.' MEG_fields{i} '(f);']);
        end
    end
    
    %Filter trials with selected Betalevel
    numoftrials = num2str(length(data.trialNum));
    f = data.stim.beta == str2double(betalevel_input); %filter to select trials for certain beta level
    
    dat_fields  = fieldnames(data);
    for i = 1:length(dat_fields)
        if eval(['length(data.' dat_fields{i} ') ==' numoftrials])
            eval(['data.' dat_fields{i} ' = data.' dat_fields{i} '(f);']);
        end
    end
    
    stim_fields = fieldnames(data.stim);
    for i = 1:length(stim_fields)
        if eval(['length(data.stim.' stim_fields{i} ') ==' numoftrials])
            eval(['data.stim.' stim_fields{i} ' = data.stim.' stim_fields{i} '(f);']);
        end
    end
    
    MEG_fields = fieldnames(data_MEG_perTD);
    for i = 1:length(MEG_fields)
        if eval(['length(data_MEG_perTD.' MEG_fields{i} ') ==' numoftrials])
            eval(['data_MEG_perTD.' MEG_fields{i} ' = data_MEG_perTD.' MEG_fields{i} '(f);']);
        end
    end
    
    %Check MEG data for NaNs (the case in S11, 300ms, 600ms) and remove
    %corresponding trials
    numoftrials = num2str(length(data.trialNum));
    NaNtrial = [];
    for i_trial = 1:length(data_MEG_perTD.trial)
        NaNfind_dataMEG{i_trial} = find(isnan(data_MEG_perTD.trial{i_trial}));
        if ~isempty(NaNfind_dataMEG{i_trial})
            NaNtrial(1,:) = i_trial;
        end
    end
    
    if ~isempty(NaNtrial)
        disp(['NaN entries found - removing corrupted trials [' num2str(NaNtrial) ']']);
        f = ~(1:length(data.trialNum) == NaNtrial); %filter to select trials for certain beta level

        dat_fields  = fieldnames(data);
        for i = 1:length(dat_fields)
            if eval(['length(data.' dat_fields{i} ') ==' numoftrials])
                eval(['data.' dat_fields{i} ' = data.' dat_fields{i} '(f);']);
            end
        end

        stim_fields = fieldnames(data.stim);
        for i = 1:length(stim_fields)
            if eval(['length(data.stim.' stim_fields{i} ') ==' numoftrials])
                eval(['data.stim.' stim_fields{i} ' = data.stim.' stim_fields{i} '(f);']);
            end
        end

        MEG_fields = fieldnames(data_MEG_perTD);
        for i = 1:length(MEG_fields)
            if eval(['length(data_MEG_perTD.' MEG_fields{i} ') ==' numoftrials])
                eval(['data_MEG_perTD.' MEG_fields{i} ' = data_MEG_perTD.' MEG_fields{i} '(f);']);
            end
        end
    end       
    
    %Store TD-wise data in proxy var and delete computational var
    data_MEG{i_tonedur} = data_MEG_perTD;
    data_MEG{i_tonedur}.behav = data;
    data_MEG{i_tonedur}.stim = data.stim;
    data_MEG_perTD =[];
    data = [];    
end

%Append data across TD
cfg = [];
cfg.keepsampleinfo = 'yes';
data_AllTD = ft_appenddata(cfg,data_MEG{1},data_MEG{2},data_MEG{3}); %appends .trial and .time
%Manually transfer and append rest
data_AllTD.hdr = data_MEG{1}.hdr;
data_AllTD.fsample = data_MEG{1}.fsample;

%Create empty subfields
behav_fields = fieldnames(data_MEG{1}.behav);
for i = 1:length(behav_fields)
        eval(['data_AllTD.behav.' behav_fields{i} ' = []']);
end
stim_fields = fieldnames(data_MEG{1}.stim);
for i = 1:length(stim_fields)
        eval(['data_AllTD.stim.' stim_fields{i} ' = []']);
end

%Fill subfields with appended TD-wise info
for i_tonedur = 1:length(tonedur_text)
    trialnum = num2str(length(data_MEG{i_tonedur}.trial));
    behav_fields = fieldnames(data_MEG{i_tonedur}.behav);
    for i = 1:length(behav_fields)   
        if eval(['length(data_MEG{i_tonedur}.behav.' behav_fields{i} ') ==' trialnum])
            eval(['data_AllTD.behav.' behav_fields{i} ...
                ' = [data_AllTD.behav.' behav_fields{i} ...
                ', data_MEG{i_tonedur}.behav.' behav_fields{i} '];']);
        end
    end
    stim_fields = fieldnames(data_MEG{i_tonedur}.stim);
    for i = 1:length(stim_fields)   
        if eval(['length(data_MEG{i_tonedur}.stim.' ...
                stim_fields{i} ') ==' trialnum])
            eval(['data_AllTD.stim.' stim_fields{i} ...
                ' = [data_AllTD.stim.' stim_fields{i} ...
                ', data_MEG{i_tonedur}.stim.' stim_fields{i} '];']);
        end
    end    
end

%% 2) Determine analysis parameters and find samples corresponding to respective tones
%Define time window parameters
nTrials  = length(data_AllTD.trial);
nSensors = size(data_AllTD.trial{1}, 1 );
samplingFreq = data_AllTD.fsample;

%TD-specific parameters, looped over TD
for i_tonedur = 1:length(tonedur_text)
    
    f_tonedur = data_AllTD.stim.toneDur == str2double(tonedur_text{i_tonedur}); %filter to select trials for certain tonedur
    index_TonedurTrials{i_tonedur} = find(f_tonedur);
    nTrials_pertonedur{i_tonedur}  = length(find(f_tonedur));
    
    toneDur_inSecs{i_tonedur}  = str2num(tonedur_text{i_tonedur});
    nSamplesPerTone{i_tonedur} = toneDur_inSecs{i_tonedur} * samplingFreq;
    
    nSamplesPerSeries = length(data_AllTD.trial{index_TonedurTrials{i_tonedur}(1)}(1,:));
    
    %Define starting point (i.e., first tone)
    t_start_ind(1) = find(data_AllTD.time{index_TonedurTrials{i_tonedur}(1)} == 0);
    t_end_ind(1)   = t_start_ind(1) + nSamplesPerTone{i_tonedur}  - 1;
    
    %Read out start and stop samples for each tone (except the first)
    for i_tone = 2:34 %loop across tones
        t_start_ind(i_tone) = t_end_ind(i_tone-1) + 1; %Start index is end index of previous tone +1
        t_end_ind(i_tone)   = t_start_ind(i_tone) + nSamplesPerTone{i_tonedur} - 1;%End index/final data point of i_tone is start+number of samples per tone -1
    end
    
    %Initialize 3D MEG arrays for storage
    MEGdata_p33_perTonedur{i_tonedur} = zeros(nSensors, nSamplesPerTone{i_tonedur}, nTrials_pertonedur{i_tonedur});
    MEGdata_p34_perTonedur{i_tonedur} = MEGdata_p33_perTonedur{i_tonedur};
    
    %copy trial wise info (i.e., all samples for selected tone for all channels and all trials into data arrays
    trialcounter = 0;
    for i_trial = index_TonedurTrials{i_tonedur} 
        trialcounter = trialcounter + 1;
        MEGdata_p33_perTonedur{i_tonedur}(:, :, trialcounter) = ...
            data_AllTD.trial{i_trial}(:, t_start_ind(toneIndex) : t_end_ind(toneIndex));        
        MEGdata_p34_perTonedur{i_tonedur}(:, :, trialcounter) = ...
            data_AllTD.trial{i_trial}(:, t_start_ind(toneIndex + 1) : t_end_ind(toneIndex + 1));   
    end
    
    clear f_tonedur t_start_ind t_end_ind
    
end

%2Combine data across TD
MEGdata_p33 = nan(nSensors, nSamplesPerTone{3}, nTrials);
MEGdata_p34 = MEGdata_p33;

for i_tonedur = 1:length(tonedur_text)    
    MEGdata_p33(:,1:size(MEGdata_p33_perTonedur{i_tonedur},2),index_TonedurTrials{i_tonedur}(1:end)) = ...
        MEGdata_p33_perTonedur{i_tonedur};   
    %Output: matrix with Chan*Samples*Trials; First 90samples/150ms are
    %filled with data for all trials
    
    MEGdata_p34(:,1:size(MEGdata_p34_perTonedur{i_tonedur},2),index_TonedurTrials{i_tonedur}(1:end)) = ...
        MEGdata_p34_perTonedur{i_tonedur};   
end

%% 3) Get predicted final tones based on each predictive_sequencerange ***
%for each trial, compute/read out the tone pitch predicted given the data so far (p*34)
for i_k = 1:length(predictive_sequencerange)
   
    series_start_ind = toneIndex - predictive_sequencerange(i_k) + 1; %beginning of tone sequence part used for prediction
    series_end_ind   = toneIndex; %end of tone sequence part used for prediction

    for i_trial = 1:nTrials
        series_predp34_discretized(i_trial) = data_AllTD.stim.logf_pred(i_trial); %log(p*34) discretized
    end
end

%% 4) Compute avg ERF activity per time window for current tone
%Define windows
win_size    = 30;
win_overlap = 0;
nSensors = size(data_AllTD.trial{1},1);

%TD-specific parameters
for i_tonedur = 1:length(tonedur_text)
    toneDur_inSecs{i_tonedur}  = str2num(tonedur_text{i_tonedur});
    nSamplesPerTone{i_tonedur} = toneDur_inSecs{i_tonedur} * samplingFreq;
    
    %Define number, start and end sample of window per tone
    windows{i_tonedur} = [1 win_size];
    while windows{i_tonedur}(end,end) < nSamplesPerTone{i_tonedur}
        windows{i_tonedur} = [windows{i_tonedur}; ...
            windows{i_tonedur}(end,:) + (win_size - win_overlap)];
    end
    
    if windows{i_tonedur}(end,end) > nSamplesPerTone{i_tonedur}
        windows{i_tonedur}(end,:) = [];
    end
    
    windows_inms{i_tonedur} = (windows{i_tonedur} / samplingFreq) * 1000;
    
end

%Compute avg ERF for selected tones (p33, p34)
for i_tonedur = 1:length(tonedur_text)    
    for i_win = 1:size(windows{3},1) %Take max window number (non-present data will remain NaN)
        
        ind_start = windows{3}(i_win, 1);
        ind_end   = windows{3}(i_win, 2);
        
        ERF_p33_win{i_win} = squeeze( nanmean( MEGdata_p33(:, ind_start:ind_end, :), 2) );
        ERF_p34_win{i_win} = squeeze( nanmean( MEGdata_p34(:, ind_start:ind_end, :), 2) );
        
        %Output: channel*trial matrix per window with 1 ERF values averaged
        %across all samples per window
    end
end

%% 5) Compute measures of association between ERF windows and prediction / prediction error at each sensor

for i_win = 1:size(windows{1},1) %Take min window number (i.e., only windows where there is data from all TD trials)
    for i_sensor = 1:nSensors
        
        %%% for prediction at tone 33 (i.e., with ERF of current tone (33)) %%%
        dv = ERF_p33_win{i_win}(i_sensor,:)';  % mean ERF window at this sensor for each trial
        iv_pred = series_predp34_discretized(:)'; %discretized predicted tone pitch for each trial
        
        % linear regression between ERF (per channel, time window) and predicted tone pitch
        stats = regstats(dv, iv_pred, 'linear', 'tstat');
        
        pred_t1{i_k}{i_win}(i_sensor) = stats.tstat.t(2);
        pred_t1_stats{i_k}{i_win}{i_sensor} = stats;
        pred_t1_p{i_k}{i_win}(i_sensor) = stats.tstat.pval(2);        

    end
end

%% 6) Save variables
mkdir([path_save]);
savefile = [path_save sub '_Pred_SeqBeta' betalevel_input '.mat'];

save(savefile, 'pred_t1', 'pred_t1_stats', 'pred_t1_p', ... 
               'predictive_sequencerange', 'betalevel_input', 'toneIndex', ...
               'ERF_p33_win', ...
               'series_predp34_discretized');
end                