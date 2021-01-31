function NASTD_MEG_Pred_CompPred_Subs...
    (sub, tonedur_text, ...
    predictive_sequencerange, toneIndex, ...
    plot_figs, pval_plotting, save_figs, ...
    paths_NASTD_MEG)

%Aim: Compute and save neural prediction effect (How is neural activity at tone 33
%modulated by the expected value of tone 34, which itself depends on the
%previous tone sequence input?)

%% 0) Specify vars, paths, and setup fieldtrip
addpath(genpath(paths_NASTD_MEG.ScriptsDir));
NASTD_MEG_SubInfo

%Tone duration condition for data load-in
if strcmp(tonedur_text,'0.15')
    tonedur_title = '0.15sTD';
elseif  strcmp(tonedur_text,'0.3')
    tonedur_title = '0.3sTD';
elseif strcmp(tonedur_text,'0.6')
    tonedur_title = '0.6sTD';
end

%Data input
loadfile_MEG = [si.path_ica  'data_clean.mat'];%path to single-subject MEG data
loadfile_behav = si.path_behav;%path to single-subject behavioral/stimulus data

%Data output
path_save = [paths_NASTD_MEG.Current_outputdata 'PredEffect/' tonedur_title '/' sub '/'];

%% 1) Load and prepare input data
%Load preprocessed MEG data
disp('loading...')
load(loadfile_MEG);
load(loadfile_behav);

%Extract clean trials
data_MEG = extract_trials(data_clean, sub, str2double(tonedur_text));

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
dat_MEG_fields  = fieldnames(data_MEG);
for i = 1:length(dat_MEG_fields)
    if eval(['length(data_MEG.' dat_MEG_fields{i} ') ==  ' numoftrials])
        eval(['data_MEG.' dat_MEG_fields{i} ' = data_MEG.' dat_MEG_fields{i} '(si.good_trials);']);
    end
end

%Filter trials with selected ToneDur
numoftrials = num2str(length(data.trialNum));

filter_tonedur = data.stim.toneDur == str2double(tonedur_text); %filter to select trials for certain tone dur

dat_fields  = fieldnames(data);
for i = 1:length(dat_fields)
    if eval(['length(data.' dat_fields{i} ') ==' numoftrials])
        eval(['data.' dat_fields{i} ' = data.' dat_fields{i} '(filter_tonedur);']);
    end
end

stim_fields = fieldnames(data.stim);
for i = 1:length(stim_fields)
    if eval(['length(data.stim.' stim_fields{i} ') ==' numoftrials])
        eval(['data.stim.' stim_fields{i} ' = data.stim.' stim_fields{i} '(filter_tonedur);']);
    end
end

MEG_fields = fieldnames(data_MEG);
for i = 1:length(MEG_fields)
    if eval(['length(data_MEG.' MEG_fields{i} ') ==' numoftrials])
        eval(['data_MEG.' MEG_fields{i} ' = data_MEG.' MEG_fields{i} '(filter_tonedur);']);
    end
end

%Append behavioral / stim data to MEG data
data_MEG.behav = data;
data_MEG.stim  = data.stim;

%% 2) Determine analysis parameters and find samples corresponding to respective tones
%Define time window parameters
nTrials  = length(data_MEG.trial);
nSensors = size( data_MEG.trial{1}, 1 );

samplingFreq    = data_MEG.fsample;
toneDur_inSecs  = str2num(tonedur_text);
nSamplesPerTone = toneDur_inSecs * samplingFreq;

nSamplesPerSeries = length(data_MEG.trial{1}(1,:));

%Define starting point (i.e., start of first tone)
t_start_ind(1) = find(data_MEG.time{1} == 0);
t_end_ind(1)   = t_start_ind(1) + nSamplesPerTone  - 1;

%Read out start and stop samples for each tone (except the first)
for i_tone = 2:34 %loop across tones
    t_start_ind(i_tone) = t_end_ind(i_tone-1) + 1; %Start index is end index of previous tone +1
    t_end_ind(i_tone)   = t_start_ind(i_tone) + nSamplesPerTone - 1;%End index/final data point of i_tone is start+number of samples per tone -1
end

%Confirm timeline
% [data_MEG.time{1}(t_start_ind); data_MEG.time{1}(t_end_ind)]'
% t_end_ind-t_start_ind

%Initialize 3D MEG arrays for storage
MEGdata_p33 = zeros(nSensors, nSamplesPerTone, nTrials);

%Copy trial wise info (i.e., all samples for selected tone for all channels and all trials into data arrays
for i_trial = 1:nTrials
    MEGdata_p33(:, :, i_trial) = ...
        data_MEG.trial{i_trial}(:, t_start_ind(toneIndex) : t_end_ind(toneIndex));
end

%% 3) Get predicted final tones based on each predictive_sequencerange ***
%for each trial, compute/read out the tone pitch predicted given the data so far (p*34)
for i_k = 1:length(predictive_sequencerange)
    series_start_ind = toneIndex - predictive_sequencerange(i_k) + 1; %beginning of tone sequence part used for prediction
    series_end_ind   = toneIndex; %end of tone sequence part used for prediction
    
    for i_trial = 1:nTrials
        series_predp34_discretized(i_trial) = data_MEG.stim.logf_pred(i_trial); %log(p*34) discretized
    end
end

%% 4) Compute avg ERF activity per time window for Toneindex tone (p33)
%Define windows
win_size    = 30;
win_overlap = 0;

samplingFreqample = data_MEG.fsample;
toneDur_inSecs  = str2num(tonedur_text);
nSamplesPerTone = toneDur_inSecs * samplingFreqample;

%Define number, start and end sample of window per tone
windows = [1 win_size];
while windows(end,end) < nSamplesPerTone
    windows = [windows; windows(end,:) + (win_size - win_overlap)];
end
if windows(end,end) > nSamplesPerTone %windows within single tone
    windows(end,:) = [];
end
windows_inms = (windows / samplingFreqample) * 1000;

nSensors = size(data_MEG.trial{1},1);

%Compute avg ERF for selected tones (p33, p34)
for i_win = 1:size(windows,1)
    ind_start = windows(i_win, 1);
    ind_end   = windows(i_win, 2);
    ERF_p33_win{i_win} = squeeze( mean( MEGdata_p33(:, ind_start:ind_end, :), 2) );
    %Output: channel*trial matrix per window with 1 ERF values averaged across all samples per window
end

%% 5) Compute measures of association between ERF windows and prediction at each sensor
for i_k = 1:length(predictive_sequencerange)
    for i_win = 1:size(windows,1)
        for i_sensor = 1:nSensors
            if predictive_sequencerange(i_k) > 1
                %only if the currently focused tone is not the first tone
                %(since for first tone we can't really predict anything)
                
                %%% for prediction at tone 33 (i.e., with ERF of current tone (33))
                dv = ERF_p33_win{i_win}(i_sensor,:)';  % mean ERF window at this sensor for each trial
                iv_pred = series_predp34_discretized(:)'; %discretized predicted tone pitch for each trial
                
                % linear regression between ERF (per channel, time window) and predicted tone pitch
                stats = regstats(dv, iv_pred, 'linear', 'tstat');
                
                pred_t1{i_k}{i_win}(i_sensor) = stats.tstat.t(2);
                pred_t1_stats{i_k}{i_win}{i_sensor} = stats;
                pred_t1_p{i_k}{i_win}(i_sensor) = stats.tstat.pval(2);
            end
        end
    end
end

%% 7) Save variables
mkdir([path_save]);
savefile = [path_save sub '_PredEffect_' tonedur_title '.mat'];

save(savefile, 'pred_t1', 'pred_t1_stats', 'pred_t1_p', ...
    'predictive_sequencerange', 'toneIndex', ...
    'ERF_p33_win', ...
    'series_predp34_discretized');

%% 8) Plot topoplot (showing uncorrected p-values per sensor for each time window
if plot_figs == 1
    %Load in label file
    load([paths_NASTD_MEG.ScriptsDir 'MEG_sensor_setup_272/label272.mat'] ); %file with CTF sensor labels for 272 sensors
    
    %Find zlimits for constant plotting across time windows
    for i_k = 1:length(predictive_sequencerange)
        
        dv_absmax = 0;
        for i_win = 1:size(windows,1)
            dv_win_absmax = max( abs( pred_t1{i_k}{i_win} ));
            if dv_win_absmax > dv_absmax
                dv_absmax = round(dv_win_absmax);
            end
        end
        
        %Prepare struct for topoplot for current window and plot
        for i_win = 1:size(windows,1)
            
            %Set up title information
            w1 = num2str( 1000 * windows(i_win, 1) / samplingFreq , 3 );
            w2 = num2str( 1000 * windows(i_win, 2) / samplingFreq , 3 );
            win_title = ['time window = [' w1 ' ms - ' w2 ' ms]'];
            
            assoc_title = 'Prediction effect (t-stat from linear regression; ';
            pval = pred_t1_p{i_k}{i_win};
            title_text = ['evaluated at tone #' num2str(toneIndex) ...
                ', prediction using k = ' num2str(predictive_sequencerange(i_k)) ' tones'];
            
            clear dat
            dat.dimord = 'chan_time';
            dat.label  = label;
            dat.time   = 0;
            
            cfg = [];
            cfg.layout    = 'CTF275.lay';
            cfg.comment   = 'no';
            cfg.colorbar  = 'yes';
            cfg.zlim      = [-dv_absmax, dv_absmax];
            
            cfg.highlight           = 'on';
            cfg.highlightcolor      = [0.9 0.9 0.9];
            cfg.highlightsymbol     = '.';
            cfg.highlightsize       = 30;
            
            cfg.highlightchannel    = find( pval <= pval_plotting);
            
            dv = pred_t1{i_k}{i_win}';
            dat.avg = dv;
            
            h = figure;
            set(gcf,'units','normalized','outerposition',[0 0 1 1])
            ft_topoplotER(cfg, dat);
            
            title({[sub ' ' assoc_title ' p < ' num2str(pval_plotting) ...
                ', uncorrected)'], ['TD: ' tonedur_title '; ERF ' win_title], title_text})
            
            %Save topoplot
            if save_figs    == 1
                path_fig = ([paths_NASTD_MEG.Current_outputfig 'PredEffect_uncorrected/' tonedur_title '/' sub '/']);
                mkdir(path_fig);
                filename     = [sub '_PredEffect_uncorrected_' tonedur_title '0_TW' w1 '-' w2 '.png'];
                figfile      = [path_fig filename];
                saveas(gcf, [figfile], 'png'); 
                delete(h);
            end
            
        end
    end
end