function NASTD_MEG_Preproc_ICAcalculation...
    (sub, paths_NASTD_MEG)

%Aim: Compute ICA on MEG raw data

%% 0) Specify vars, paths, and setup fieldtrip
addpath(genpath(paths_NASTD_MEG.ScriptsDir));
addpath(paths_NASTD_MEG.FieldTrip);
ft_defaults

NASTD_MEG_SubInfo

mkdir(si.path_data_sub, 'ICA');
path_outputdata = [si.path_data_sub 'ICA/'];

savefile = [path_outputdata 'data_ica.mat']; %subject specific output file 

%% 1) load events

load([si.path_events 'events.mat']);


%% 2) load and filter block-length meg data segments


for i_block = 1:length(si.blocks)

    cfg = [];
    cfg.dataset = si.path_raw{i_block};
    cfg.channel = {'MEG'};
    d = ft_preprocessing(cfg);    

%% 3) Remove excessive data at the end of the block 
    
%     % take out trailing zeros    
%     [firstEmptySample d] = find_empty_MEG_data(d);

    % chop off recording after 0 sec following final button press 
    nSensors = size(d.trial{1}, 1);
    
    fs = d.fsample;
    postResponseDur_inSecs = 0;
    postResponseDur_inSamples = postResponseDur_inSecs * fs;

%     if si.isConf
%         last_resp_sample = event_resp3{i_block}(end).sample;
%     else
%         last_resp_sample = event_resp2{i_block}(end).sample;
%     end
    
    last_resp_sample = event_resp2{i_block}(end).sample;

    last_sample = last_resp_sample + postResponseDur_inSamples;
    newTrial = zeros(nSensors, last_sample);
    for i_sensor = 1:nSensors
        newTrial(i_sensor, :) = d.trial{1}(i_sensor, 1:last_sample);
    end
    
    d.trial{1} = newTrial;
    d.time{1}  = d.time{1}(1:last_sample);    

    d.sampleinfo = [1, length(d.time{1})];

    
%% 4) Apply preprocessing filters for ICA

    cfg = [];
    cfg.continuous = 'yes';
    cfg.demean     = 'yes';
    cfg.detrend    = 'yes';

    cfg.bsfilter   = 'yes';
    cfg.bsfreq     = [58 62; 118 122];
    cfg.bsfiltord  = 4;
    cfg.bsfilttype = 'but';

    cfg.hpfilter   = 'yes';
    cfg.hpfreq     = .05;
    cfg.hpfiltord  = 3;
    cfg.hpfilttype = 'but';
    
    data{i_block} = ft_preprocessing(cfg, d);

end


%% 5) Concatenate blocks, treating them as trials

data_all.hdr     = data{1}.hdr;
data_all.label   = data{1}.label;
data_all.fsample = data{1}.fsample;
data_all.grad    = data{1}.grad;
data_all.cfg     = data{1}.cfg;

data_all.time       = [];
data_all.trial      = [];
data_all.sampleinfo = [];

for i_block = 1:length(si.blocks)
    data_all.time       = [data_all.time  data{i_block}.time];
    data_all.trial      = [data_all.trial data{i_block}.trial];
    
    if i_block == 1
        data_all.sampleinfo = [data_all.sampleinfo;  data{i_block}.sampleinfo];
    else
        data_all.sampleinfo = [data_all.sampleinfo;  data{i_block}.sampleinfo ...
            + data_all.sampleinfo(end,end)];
    end
end

%% 6) Calculate ICA
cfg        = [];
cfg.method = 'runica'; % try also pca, with cfg.runica.pca=32;
% cfg.runica.pca = 32;
% cfg.runica.pca = 16;
data_ica = ft_componentanalysis(cfg, data_all);

%% 7) Save the preprocessed data
% save(savefile, 'data_ica', 'cfg', 'tf', 'tf_hrs', '-v7.3');

end