function NASTD_MEG_Preproc_ICAremove...
    (sub, paths_NASTD_MEG)

%Aim: Remove ICA-components containing artifacts and save cleaned data

%% 0) Specify vars, paths, and setup fieldtrip
addpath(genpath(paths_NASTD_MEG.ScriptsDir));
addpath(paths_NASTD_MEG.FieldTrip);
ft_defaults

NASTD_MEG_SubInfo

loadfile_ica = [si.path_ica  'data_ica.mat'];
savefile = [si.path_ica 'data_clean.mat'];

path_fig = [si.path_data_sub 'Figs/'];

%% 1) Load ICA and MEG data

% load data_ica.mat
disp('loading...')
load(loadfile_ica);
disp('done.')

% load selected ICA components for removal
run([si.path_ica 'selected_artifacts.m']);
artifact_components = ac_all;

% load events
load([si.path_events 'events.mat']);

%% 2) load and filter block-length MEG data segments
for i_block = 1:length(si.blocks)
    
    cfg = [];
    cfg.dataset = si.path_raw{i_block};
    cfg.channel = {'MEG'};
    d = ft_preprocessing(cfg);
    
    
    % remove excessive data at the end of the block
    
    %     % take out trailing zeros
    %     [firstEmptySample d] = find_empty_MEG_data(d);
    
    % chop off recording after 2 sec following final button press
    nSensors = size(d.trial{1}, 1);
    
    fs = d.fsample;
    postResponseDur_inSecs = 2;
    postResponseDur_inSamples = postResponseDur_inSecs * fs;
    
    last_resp_sample = event_resp2{i_block}(end).sample;
    
    last_sample = last_resp_sample + postResponseDur_inSamples;
    
    % for subject S11, the last few trials of the last block were omitted accidentally
    % so for this subject, we define the last sample manually to ensure
    % maximum retention of data
    if strcmp(sub,'S11') && i_block == 12
        last_sample = 235978;
    end
    
    newTrial = zeros(nSensors, last_sample);
    for i_sensor = 1:nSensors
        newTrial(i_sensor, :) = d.trial{1}(i_sensor, 1:last_sample);
    end
    
    d.trial{1} = newTrial;
    d.time{1}  = d.time{1}(1:last_sample);
    
    d.sampleinfo = [1, length(d.time{1})];
    
    
    %% 3) apply filters for ICA    
    cfg = [];
    cfg.continuous = 'yes';
    cfg.demean     = 'yes';
    cfg.detrend    = 'yes';
    
    cfg.bsfilter   = 'yes';
    cfg.bsfreq     = [58 62; 118 122];
    cfg.bsfiltord  = 4;
    cfg.bsfilttype = 'but';
    
    % do not apply HP filter for cleaned data
    cfg.hpfilter   = 'no';
    %     cfg.hpfreq     = .05;
    %     cfg.hpfiltord  = 3;
    %     cfg.hpfilttype = 'but';
    
    data{i_block} = ft_preprocessing(cfg, d);
    
end

%% 4) concatenate blocks, treating them as trials

data_all.hdr     = data{1}.hdr;
data_all.label   = data{1}.label;
data_all.fsample = data{1}.fsample;
data_all.grad    = data{1}.grad;
data_all.cfg     = data{1}.cfg;

data_all.time       = [];
data_all.trial      = [];
data_all.sampleinfo = [];

for i_block = 1:length(si.blocks)
    data_all.time       = [data_all.time         data{i_block}.time];
    data_all.trial      = [data_all.trial        data{i_block}.trial];
    
    if i_block == 1
        data_all.sampleinfo = [data_all.sampleinfo;  data{i_block}.sampleinfo];
    else
        data_all.sampleinfo = [data_all.sampleinfo;  data{i_block}.sampleinfo + data_all.sampleinfo(end,end)];
    end
end

%% 5) Apply the ICA correction

cfg           = [];
cfg.component = artifact_components; % to be removed component(s)

data_clean    = ft_rejectcomponent(cfg, data_ica, data_all);

%% 6) Manually review the ICA correction
% this section of code is meant to be run manually as a follow-up to manual
% inspection of the ICA data for artifact selection

screen_size = get(0,'ScreenSize');
screen_l = screen_size(3);
screen_h = screen_size(4);

% butterfly plot of original, preprocessed data
cfg            = [];
cfg.viewmode   = 'butterfly';
cfg.continuous = 'no';
% cfg.channel    = {'MEG'};

ft_databrowser(cfg, data_all);
set(gcf, 'Position', [1, 1, (screen_l/2)-1, screen_h]);

% butterfly plot of new, rejected artifact data
cfg            = [];
cfg.viewmode   = 'butterfly';
cfg.continuous = 'no';

ft_databrowser(cfg, data_clean);
set(gcf, 'Position', [(screen_l/2), 1, screen_l/2, screen_h]);

return

%% 7) save the ICA correction and illustrating figures

% disp('saving...')
% save(savefile, 'data_clean', '-v7.3');
% save([si.path_ica 'artifact_components.mat'], 'artifact_components');
% disp('done.')

for i_block = 1:length(si.blocks)
    h = figd(15,1); hold on;
    plot(data_all.time{i_block}, data_all.trial{i_block}')
    
    xlabel('time (sec)')
    ylabel('T')
    title(['pre-processed data, block ' si.blocks{i_block}])
    
    format = 'png';
    figfile = [path_fig 'block_' si.blocks{i_block} '_pre.' format];
    save_fig(h, figfile, 1, 2, format);
    delete(h);
end

for i_block = 1:length(si.blocks) %1:length(si.blocks)
    h = figd(15,1); hold on;
    plot(data_clean.time{i_block}, data_clean.trial{i_block}')
    
    xlabel('time (sec)')
    ylabel('T')
    title(['ICA cleaned data, block ' si.blocks{i_block}])
        
    format = 'png';
    figfile = [path_fig 'block_' si.blocks{i_block} '_with_ICA.' format];
    save_fig(h, figfile, 1, 2, format);
    delete(h);
end

end
