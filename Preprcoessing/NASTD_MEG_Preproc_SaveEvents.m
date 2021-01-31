function NASTD_MEG_Preproc_SaveEvents...
    (sub, paths_NASTD_MEG)

%Aim: Read out events from raw MEG data, categorize them, adjust them to
%presentation delay, and then save events

%% 0) Specify vars, paths, and setup fieldtrip
addpath(genpath(paths_NASTD_MEG.ScriptsDir));
addpath(paths_NASTD_MEG.FieldTrip);
ft_defaults

NASTD_MEG_SubInfo

path_outputdata = [si.path_data_sub 'events/'];
mkdir(path_outputdata)
savefile = 'events.mat';

%% 1) Read in events per block
for i_block = 1:length(si.blocks)

    % read the triggers from the raw data
    event = ft_read_event(si.path_raw{i_block});

    % since the first trigger of 255 is read in as [], fix it...
    event(1).value = 255;

    for i = 1:length(event)
        if isempty(event(i).value), event(i).value = Inf; end
    end
    
    event_all{i_block} = event;
    
    % filter triggers by trigger ID
    % 11:15 --> onset of auditory stimulus with betaID = 1:5; 
    % 16    --> onset of REST trial 
    event_stim{i_block} = ft_filter_event(event, 'type', 'UPPT002', 'value', 11:16);
    event_stim_orig{i_block} = event_stim{i_block};
    
    % filter triggers by response ID
    % 21:25 / 29 --> response 1 (final tone prob) / failure to enter response
    event_resp1{i_block} = ft_filter_event(event, 'type', 'UPPT002', 'value', 20:29);
    
    % 31:35 / 39 --> response 2 (trend strength)  / failure to enter response
    event_resp2{i_block} = ft_filter_event(event, 'type', 'UPPT002', 'value', 30:39);

    % 41:45 / 49 --> response 3 (confidence)      / failure to enter response    
    event_resp3{i_block} = ft_filter_event(event, 'type', 'UPPT002', 'value', 40:49);
   

    %% Update stim onsets
    load([si.path_events 'stim_onset_adj_' num2str(i_block) '.mat']);

    for i = 1:length(event_stim{i_block})

        if event_stim{i_block}(i).sample ~= original_event_onsets(i)
            disp('old trigger values don''t match up')
            keyboard
            
        else
            event_stim{i_block}(i).sample = new_event_onsets(i);
            
        end

    end
    
end

%% save data
% save([path_outputdata, savefile], ...
%     'event_all', 'event_stim', 'event_stim_orig', 'event_resp1', 'event_resp2', 'event_resp3');
