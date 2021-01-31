function NASTD_MEG_ERF_CompM100SpatialFilter...
    (i_tonedur, ...
    param, paths_NASTD_MEG)
%Aim: Compute individual spatial filter (sensor level) for M100
%1) Compute individual M100 ERF
%2) Compute linear weights for each sensor within M100
%(i.e., how much does each respective sensor contribute to M100 signal)?
%3) Apply spatial filter to trial-wise data
%4) plot summary figure showing single-subject and group-level filter
%weights

%% 1. Compute individual ERF for each TD and place it in joint FT-struct
path_fig = [paths_NASTD_MEG.Current_outputfig 'M100SpatialFilter/Ssub/' param.tonedur_text{i_tonedur} 'sTD/'];
path_output = [paths_NASTD_MEG.Current_outputdata 'Data_M100SpatialFilter/Ssub/' param.tonedur_text{i_tonedur} 'sTD/'];
mkdir([path_output]);
mkdir(path_fig)

IncludedTones = param.M100.IncludedTones;

sub_counter = 0;
for i_sub = 1:length(param.subs) %All subs
    
    sub_counter = sub_counter+1;
    disp([param.subs{sub_counter} ' - TD: ' param.tonedur_text{i_tonedur} ' - Start M100 spatial filter computation'])
    
    %Cut out trials per subject for specific TD, read out trial indices for
    %specific subconditions, and provide a file ID for identification.
    MEGactivity_perTrial = ...
        NASTD_MEG_ERF_DefineTrialsperTD(...
        param.subs{i_sub}, ...
        param.tonedur_text{i_tonedur}, ...
        param.samplefreq, ...
        paths_NASTD_MEG);
    
    %Output: Struct with cleaned MEG data for all trials for specific TD, including stim/trial info
    MEGactivity_perTrial.trialindices.predp34_low = ...
        find(MEGactivity_perTrial.stim.predID == -1); %p*34 = low
    MEGactivity_perTrial.trialindices.predp34_high = ...
        find(MEGactivity_perTrial.stim.predID == 1); %p*34 = high
    
    data_length_sec = length(MEGactivity_perTrial.trial{1})/param.samplefreq;
    MEGactivity_perTrial.ID = ...
        [param.subs{i_sub} '_' ...
        param.tonedur_text{i_tonedur} 'sTD_' ...
        num2str(data_length_sec) 'sDataLength'];
    
    %Determine tone-dependent begin/end time
    nSamplesPerTone = str2double(param.tonedur_text{i_tonedur}) * param.samplefreq;
    index_T1start = find(MEGactivity_perTrial.time{1} == 0); %determine start of first tone
    Samples_AllToneStart = [index_T1start:nSamplesPerTone:(index_T1start+33*nSamplesPerTone)];
    
    %Check for trials ocntaining NaNs (e.g., sub 6)
    allNaNtrials = [];
    for i_trial = 1:length(MEGactivity_perTrial.trial)
        if sum(isnan(MEGactivity_perTrial.trial{i_trial})) > 0
            disp(['Nan in trial: ' num2str(i_trial) ' - removing trial'])
            allNaNtrials = [allNaNtrials; i_trial];
            
        end
    end
    if ~isempty(allNaNtrials)
        alltrials = 1:length(MEGactivity_perTrial.trial);
        allnoNaNtrials = alltrials(find(alltrials ~= allNaNtrials));
        cfg = [];
        cfg.trials =  allnoNaNtrials;
        MEGactivity_perTrial = ft_selectdata(cfg,MEGactivity_perTrial);
    end
    
    %Preprocess each trial (Detrend, Demean, baseline-correct (prestim -0.5-0), and LP filter)
    cfg                     = [];
    cfg.continuous          = 'no';
    cfg.demean              = 'yes';
    cfg.detrend             = 'yes';
    cfg.baselinewindow      = [MEGactivity_perTrial.time{1}(1) MEGactivity_perTrial.time{1}(index_T1start-1)];
    cfg.lpfilter            = 'yes';% apply lowpass filter
    cfg.lpfreq              = 35; % lowpass at 35 Hz.
    ERF_perTrial_BC{i_sub}    = ft_preprocessing(cfg,MEGactivity_perTrial);
    
    %copy subfields
    ERF_perTrial_BC{i_sub}.behav = MEGactivity_perTrial.behav;
    ERF_perTrial_BC{i_sub}.stim = MEGactivity_perTrial.stim;
    ERF_perTrial_BC{i_sub}.trialindices = MEGactivity_perTrial.trialindices;
    
    clear MEGactivity_perTrial
    
    %For each trial, separate whole trial into each tone
    for i_tone = IncludedTones
        cfg             = [];
        cfg.latency     = ...
            [ERF_perTrial_BC{i_sub}.time{1}(Samples_AllToneStart(i_tone)) ...
            ERF_perTrial_BC{i_sub}.time{1}(Samples_AllToneStart(i_tone)+nSamplesPerTone-1)];
        cfg.feedback    = 'no';
        ERF_perTrial_perTone{i_tone} = ft_selectdata(cfg, ERF_perTrial_BC{i_sub});
    end
    
    %Average ERF across trials for selected tones
    time_tone1 = ERF_perTrial_perTone{1}.time{1}; %copy timeline
    for i_tone = IncludedTones
        cfg             = [];
        cfg.keeptrials  = 'no';
        ERF_AvgTrials_AvgTone{i_tone} = ft_timelockanalysis(cfg,ERF_perTrial_perTone{i_tone});
        ERF_AvgTrials_AvgTone{i_tone}.time = time_tone1; %substitute identical timeline for each tone (0:end TD)
    end
    
    %Average each sample across tones (makes only sense if multiple
    %tones selected)
    cfg = [];
    cfg.parameter       = 'avg';
    cfg.keepindividual  = 'no';
    ERF_perTrial_AvgTone = ft_timelockgrandaverage(cfg,ERF_AvgTrials_AvgTone{:});
    
    %Transform raw amplitude to power by squaring signal at each sample
    for i_sample = 1:size(ERF_perTrial_AvgTone.avg,2)
        ERF_perTrial_AvgTone.avg_squared(:,i_sample) = ...
            ERF_perTrial_AvgTone.avg(:,i_sample).^2;
    end
    
    %Place sub-field into common cell
    ERFperTone.Ssub.avg_squared{sub_counter} = ERF_perTrial_AvgTone;
    
    %     %Plot topoplot for selected M100 time window in single plot
    %     cfg                     = [];
    %     cfg.xlim                = param.M100.TW;
    %     cfg.layout              = 'CTF275';
    %     cfg.parameter           = 'avg_squared';
    %     cfg.colorbar            = 'yes';
    %     figure;
    %     title([param.subs(sub_counter)])
    %     ft_topoplotER(cfg, ERFperTone.Ssub.avg_squared{i_sub})
    
    %% 2 Compute & Apply spatial filter
    %For M100 TW (averaged across included samples), compute relative
    %sensor-contribution to total signal strength
    
    %Determine M100TW related samples
    sample_M100TW_start = find(ismember(ERF_perTrial_AvgTone.time,param.M100.TW(1)));
    sample_M100TW_stop = find(ismember(ERF_perTrial_AvgTone.time,param.M100.TW(2)));
    disp(['Selected time window: ' ...
        num2str(ERF_perTrial_AvgTone.time(sample_M100TW_start)) '-' ...
        num2str(ERF_perTrial_AvgTone.time(sample_M100TW_stop)) 's'])
    
    %Average activity across M100TW for each channel
    Avg_Squared_M100activity = mean(ERF_perTrial_AvgTone.avg_squared...
        (:,sample_M100TW_start:sample_M100TW_stop),2);
    
    %Compute relative contribution of each sensor to M100 activity
    FilterWeights_M100activity = Avg_Squared_M100activity./...
        sum(Avg_Squared_M100activity);
    
    %Apply spatial filter to single-trial data
    ERF_M100filtweight_perTrial = ERF_perTrial_BC{i_sub}; %copy subfields
    %Multiply weights with time series data for each sample and overwrite
    %trial subfield
    for i_trial = 1:size(ERF_perTrial_BC{i_sub}.trial,2)
        for i_sample = 1:size(ERF_perTrial_BC{i_sub}.time{i_trial},2)
            ERF_M100filtweight_perTrial.trial{i_trial}(:,i_sample) = ...
                ERF_perTrial_BC{i_sub}.trial{i_trial}(:,i_sample).*...
                FilterWeights_M100activity;
        end
    end
    
    %Save M100 filter weights
    savefile = [path_output  param.subs{i_sub} '_M100SpatFilt' ...
        '_FilterWeights_' param.tonedur_text{i_tonedur} 'sTD.mat'];
    save(savefile, 'FilterWeights_M100activity', '-v7.3');
    
    %Save weighted single trials
    savefile = [path_output  param.subs{i_sub} '_M100SpatFilt' ...
        '_WeightedData_' param.tonedur_text{i_tonedur} 'sTD.mat'];
    save(savefile, 'ERF_M100filtweight_perTrial', '-v7.3');    
    
    %     %Optional: plot M100 summary figure to determine subjects showing
    %     %valid M100 response
    %     %Average time domain data across trials (for plotting)
    %     cfg             = [];
    %     cfg.keeptrials  = 'no';
    %     cfg.parameter   = 'trial';
    %     TimeSeriesData_SpatFiltWeight = ft_timelockanalysis(cfg, MEGactivity_perTrial_BC{i_sub});
    %
    %     %Plot sensor weights (line and topo plot) and filtered data time course
    %     Threshold_SensDef = 0.75; %Arbitrary, just to visualize the max sensors
    %
    %     figure;
    %     set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    %
    %     s1 = subplot(2,2,1);
    %     plot(1:length(FilterWeights_M100activity),FilterWeights_M100activity,'b-','LineWidth',2)
    %     s1.XLim = [1 length(FilterWeights_M100activity)];
    %     s1.XLabel.String = 'Sensor Number';
    %     s1.YLabel.String = 'Relative Sensor Contribution';
    %     MaxSens = find(FilterWeights_M100activity > ...
    %         (max(FilterWeights_M100activity)*Threshold_SensDef));
    %     hold on;
    %     plot(MaxSens,FilterWeights_M100activity(MaxSens),'k.','MarkerSize',20)
    %     hold on;
    %     plot(s1.XLim,[max(FilterWeights_M100activity)*Threshold_SensDef ...
    %         max(FilterWeights_M100activity)*Threshold_SensDef],'k--')
    %     title({['Filter Weight per normal channel order (Visualization Thresh = ' ...
    %         num2str(Threshold_SensDef) ')']})
    %
    %     %Plot filter weight as function of ordered sensor number
    %     s2 = subplot(2,2,2);
    %     [a,b] = sort(FilterWeights_M100activity);
    %     plot(1:length(FilterWeights_M100activity),sort(FilterWeights_M100activity),'b-','LineWidth',2)
    %     title('Filter Weight per sorted channel order')
    %     s2.XLim = [1 length(FilterWeights_M100activity)+5];
    %     s2.XTickLabel = b;
    %     s2.XTick = [1:15:length(FilterWeights_M100activity)];
    %     s2.XLabel.String = 'Sensor Number (Ordered)';
    %     s2.YLabel.String = 'Relative Sensor Contribution';
    %
    %     %Plot topoplot showing filter weight for M100 window
    %     s3 = subplot(2,2,3);
    %     cfg                     = [];
    %     cfg.xlim                = param.M100.TW;
    %     cfg.layout              = 'CTF275';
    %     cfg.parameter           = 'avg_squared';
    %     cfg.colorbar            = 'yes';
    %     cfg.comment             = 'no';
    %     cfg.marker              = 'on';
    %     cfg.markersymbol        = '.';
    %     cfg.markersize          = 6;
    %     cfg.markercolor         = [0.5 0.5 0.5];
    %     cfg.highlight           = 'numbers';
    %     cfg.highlightfontsize       = 12;
    %     cfg.highlightchannel    = find(FilterWeights_M100activity > ...
    %         (max(FilterWeights_M100activity)*Threshold_SensDef));
    %     ft_topoplotER(cfg, ERFperTone.Ssub.avg_squared{i_sub})
    %     title('Squared avg amplitude for M100 TW')
    %
    %     sgtitle({[param.subs{i_sub} ' - Filter weights computed for - M100 (' ...
    %         num2str(cfg.xlim(1)) '-' num2str(cfg.xlim(2)) ...
    %         's) ']})
    %
    %     s4 = subplot(2,2,4); %Option for all GAvg all trials - no specific conditions
    %     YSize = max(mean(TimeSeriesData_SpatFiltWeight.avg));
    %     plot(1:length(TimeSeriesData_SpatFiltWeight.avg),mean(TimeSeriesData_SpatFiltWeight.avg),'r-','LineWidth',2)
    %     for i_sample = 1:length(Samples_AllToneStart)
    %         hold on;
    %         plot([Samples_AllToneStart(i_sample) Samples_AllToneStart(i_sample)] , ...
    %             [YSize*-1 YSize],'k-','LineWidth',0.5)
    %     end
    %     hold on;
    %     plot([1 length(TimeSeriesData_SpatFiltWeight.avg)] , ...
    %         [0 0],'k-','LineWidth',0.5)
    %     s4.XLabel.String = 'Time [s]';
    %     s4.XTick = 1:150:length(TimeSeriesData_SpatFiltWeight.avg);
    %     s4.XTickLabel = TimeSeriesData_SpatFiltWeight.time(1:150:end);
    %     s4.XTickLabelRotation = 270;
    %     s4.YLabel.String = 'Amplitude (weighted avg across sensors)';
    %     title('Weighted amplitude averaged across sensors')
    
    
    %Cleanup
    clear ERF_perTrial_perTone ERF_perTrial_AvgTone ERF_AvgTrials_AvgTone ...
        ERF_perTrial_AvgTone
end

%Plot Group-average topoplot of sensor-weights
for i_sub = param.M100.Subs{i_tonedur} %Only subs with valid M100 response
    
    %Load filter weights
    savefile = [path_output  param.subs{i_sub} ...
        '_M100SpatFilt_FilterWeights_' param.tonedur_text{i_tonedur} 'sTD.mat'];
    load(savefile) %FilterWeights_M100activity
    
    FilterWeights_M100activity_Allsub(:,i_sub) = FilterWeights_M100activity;
    clear FilterWeights_M100activity
end

%Set up topoplot struct
load([paths_NASTD_MEG.ScriptsDir 'MEG_sensor_setup_272/label272.mat']);
%file with CTF sensor labels for 272 sensors, called 'label'
dat.dimord = 'chan_time';
dat.label  = label;
dat.time   = 0;
dat.avg    = mean(FilterWeights_M100activity_Allsub,2); %Place GAvg sensor weights in struct

f3 = figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen

cfg                     = [];
cfg.xlim                = param.M100.TW;
cfg.layout              = 'CTF275';
cfg.parameter           = 'avg';
cfg.colorbar            = 'yes';
cfg.comment             = 'no';
cfg.zlim                = [-0.025 0.025];

title(['GAvg M100 spatial filter weights; ' param.tonedur_text{i_tonedur} ' sTD'])
ft_topoplotER(cfg, dat);

if param.plot.save
    Figtitle = (['GroupLevel_Topo_n' num2str(length(param.M100.Subs{i_tonedur})) ...
        '_' param.tonedur_text{i_tonedur} ...
        'sTD_SpatialFilterWeightsM100.png']);
    saveas(gcf, [path_fig Figtitle], 'png'); %save png version    
    close
end

end