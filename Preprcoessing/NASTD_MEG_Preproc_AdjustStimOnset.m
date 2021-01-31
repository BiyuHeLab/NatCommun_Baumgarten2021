function NASTD_MEG_Preproc_AdjustStimOnset...
    (sub, paths_NASTD_MEG)

%Aim: Read out events from raw MEG data and adjust triggers representing
%stimulus presentation to presentation delay

%% 0) Specify vars, paths, and setup fieldtrip
addpath(genpath(paths_NASTD_MEG.ScriptsDir));
addpath(paths_NASTD_MEG.FieldTrip);
ft_defaults

NASTD_MEG_SubInfo

path_outputdata = [si.path_data_sub 'events/'];
mkdir(path_outputdata)

%% 1) Read in events per block
for i_block = 1:length(si.blocks)

    path_inputdata  = si.path_raw{i_block};
    savefile        = ['stim_onset_adj_' num2str(i_block) '.mat'];

    %Read in events from raw data
    event = ft_read_event(path_inputdata);

    %since the first trigger of 255 is read in as [], fix it...
    event(1).value = 255;

    for i = 1:length(event)
        if isempty(event(i).value) 
            event(i).value = Inf; 
        end
    end

    %filter triggers by trigger ID
    %11:15 --> onset of auditory stimulus with betaID = 1:5; 
    %16    --> onset of REST trial 
    event_stim = ft_filter_event(event, 'type', 'UPPT002', 'value', 11:16);

%% 2) Read out sample values from stimulus channel and adjust them to 
%reflect the time point where the stimulus was really presented

    cfg             = [];
    cfg.dataset    = path_inputdata;
    fg.channel     = {'MEG'};
    cfg.continuous = 'yes';
    cfg.demean     = 'no';
    cfg.detrend    = 'no';

    data = ft_preprocessing(cfg);

    uadc_ind = find(strcmp(data.label, 'UADC004')); %Stimulus channel

    %Plot overview of stimulus signal
    figure;
    plot(data.trial{1}(uadc_ind,:));
    hold on;

    clear old_event new_event
    
    for i = 1:length(event_stim)

        yl = ylim;

        x = event_stim(i).sample;
        y = data.trial{1}(uadc_ind,:);

        plot([x x], yl, 'g-') %Old stimulus start sample (green)

        old_event(i) = x;

        flag_zero = 1;
        flag_neg  = 0;
        flag_pos  = 0;
        done = 0;

        % original values for triggers don't account for delays b/t sending of
        % trigger event and actual delivery of the tone sequence.
        %
        % so they briefly precede the "true" onset of auditory stimulus
        % delivery as marked by the UADC channel.
        %
        % so we want to find the sample in the data corresponding to the onset
        % of input to the UADC channel.
        %
        % so we start from the old trigger value and keep incrementing it until
        % we detect a peak that surpasses some threshold value. This we
        % consider the sample corresponding to the true onset of the auditory
        % stimulus.
        
        x = x-1;
        thresh = 10^-4;
        while ~done
            x = x+1;
            diff = y(x+1) - y(x);

            if flag_zero
                if sign(diff) < 0 && abs(diff) > thresh
                    flag_neg  = 1;
                    flag_zero = 0;
                elseif sign(diff) > 0 && abs(diff) > thresh
                    flag_pos  = 1;
                    flag_zero = 0;
                end

            elseif flag_neg
                if sign(diff) > 0 && abs(diff) > thresh
                    done = 1;
                end

            elseif flag_pos
                if sign(diff) < 0 && abs(diff) > thresh
                    done = 1;
                end
            end

        end

        new_event(i) = x;

        plot([x x], yl, 'r-') %New/Adjusted stimulus start sample (red)

    end
    
    new_event - old_event
    mean(new_event - old_event)
    max(new_event - old_event)
    min(new_event - old_event)

    title(['min = ' num2str(min(new_event - old_event)) ', max = ' num2str(max(new_event - old_event))])

    %Add 10 ms (= .01 sec * 600 samples/sec = 6 samples) to account for
    %stimulus presentation delay
    original_event_onsets = old_event;
    new_event_onsets      = new_event + 6;

%     save([path_outputdata, savefile], 'original_event_onsets', 'new_event_onsets');

end