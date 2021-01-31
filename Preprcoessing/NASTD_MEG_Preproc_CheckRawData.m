function NASTD_MEG_Preproc_CheckRawData...
    (sub, paths_NASTD_MEG)

%Aim: Visualize raw MEG data to check for bad sensors

%% 0) Specify vars, paths, and setup fieldtrip
addpath(genpath(paths_NASTD_MEG.ScriptsDir));
addpath(paths_NASTD_MEG.FieldTrip);
ft_defaults

NASTD_MEG_SubInfo

path_fig = [si.path_data_sub 'Figs/'];
mkdir(path_fig)

path_outputdata = [si.path_data_sub 'raw/'];
mkdir(path_outputdata)

savefile = [path_outputdata 'sensor_std.mat'];

%% 1) Load events and block-wise MEG data

load([si.path_events 'events.mat']);

all_trials = [];

for i_block = 1:length(si.blocks)
    
    cfg         = [];
    cfg.dataset = si.path_raw{i_block};
    cfg.channel = {'MEG'};
    d           = ft_preprocessing(cfg);
    
    
    %% 2) Remove excess empty data at the end of the block
    
    % chop off recording after 2 sec following final button press
    nSensors = size(d.trial{1}, 1);
    
    fs = d.fsample;
    postResponseDur_inSecs = 2;
    postResponseDur_inSamples = postResponseDur_inSecs * fs;
    
    last_resp_sample = event_resp2{i_block}(end).sample;
    
    last_sample = last_resp_sample + postResponseDur_inSamples;
    newTrial = zeros(nSensors, last_sample);
    for i_sensor = 1:nSensors
        newTrial(i_sensor, :) = d.trial{1}(i_sensor, 1:last_sample);
    end
    
    d.trial{1} = newTrial;
    d.time{1}  = d.time{1}(1:last_sample);
    
    % baseline correct each sensor using its mean value in first 300 samples
    for i=1:272
        d.trial{1}(i,:) = d.trial{1}(i,:) - mean( d.trial{1}(i,1:300) );
    end
    
    sensor_std_block(:,i_block) = std(d.trial{1},0,2);
    
    all_trials = [all_trials d.trial{1}];
    
    
    %     %Optionally: Preprocess data
    %     cfg = [];
    %
    %     % demean and detrend the raw recordings
    %     cfg.demean     = 'yes';
    %     cfg.detrend    = 'yes';
    %
    %
    %     % apply a filter to the data so that we take out power line noise at the 60
    %     % Hz frequency, plus its harmonic at 120 Hz
    %     cfg.bsfilter = 'yes';
    %     cfg.bsfreq   = [58 62; 118 122];
    %
    %     cfg.hpfilter  = 'yes';
    %     cfg.hpfreq    = .05;
    %     cfg.hpfiltord = 3;
    %
    %     d = ft_preprocessing(cfg, d);
    
    
    %% 3) Plot summary
    
    h = figd(15,1); hold on;
    plot(d.time{1}, d.trial{1}')
    t_end = d.time{1}(last_resp_sample);
    plot([t_end, t_end], ylim, 'k--', 'LineWidth', 2);
    xlabel('time (sec)')
    ylabel('T')
    
    format = 'png';
    figfile = [path_fig 'block_' si.blocks{i_block} '.' format];
    
    if strcmp(sub, 'S3') && i_block <= 5
        figfile = [path_fig 'block_a' si.blocks{i_block} '.' format];
    elseif strcmp(sub, 'S3') && i_block > 5
        figfile = [path_fig 'block_b' si.blocks{i_block} '.' format];
    end
    
%     save_fig(h, figfile, 1, 2, format);
%     delete(h);
end


%% 4) Save STD of each sensor across all blocks

sensor_std_all = std(all_trials,0,2);
% save(savefile, 'sensor_std_block', 'sensor_std_all');

label = d.label;
% save([path_outputdata 'label.mat'], 'label');