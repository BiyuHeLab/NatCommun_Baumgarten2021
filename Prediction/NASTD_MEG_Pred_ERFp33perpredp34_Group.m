function NASTD_MEG_Pred_ERFp33perpredp34_Group...
    (subs,  tonedur_text,...
    toneIndex, ...
    plot_figs, save_figs,...
    paths_NASTD_MEG)
%Aim: Plot timelocked ERF during p33 in predictive processing clusters 
%as a function of p*34

%% 0) Specify vars, paths, and setup fieldtrip
%Tone duration condition for data load-in
if strcmp(tonedur_text,'0.15')
    tonedur_title = '0.15sTD';
elseif  strcmp(tonedur_text,'0.3')
    tonedur_title = '0.3sTD';
elseif strcmp(tonedur_text,'0.6')
    tonedur_title = '0.6sTD';
end

%Data input
path_input = [paths_NASTD_MEG.Current_outputdata 'PredEffect/' tonedur_title '/'];
%Data output
path_output = [paths_NASTD_MEG.Current_outputdata 'ERFp33/' tonedur_title '/'];

dataMEG_p33_allsub = struct;

%% 1) Load in ERF data per sub
for i_sub = 1:length(subs)
    sub = subs{i_sub};
    NASTD_MEG_SubInfo
    %Load and prepare input data
    loadfile_MEG = [si.path_ica  'data_clean.mat'];
    loadfile_behav = si.path_behav;
    
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
    
    f = data.stim.toneDur == str2double(tonedur_text); 
    
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
    
    MEG_fields = fieldnames(data_MEG);
    for i = 1:length(MEG_fields)
        if eval(['length(data_MEG.' MEG_fields{i} ') ==' numoftrials]) 
            eval(['data_MEG.' MEG_fields{i} ' = data_MEG.' MEG_fields{i} '(f);']);
        end
    end
    
    %Append behavioral / stim data to MEG data
    data_MEG.behav = data;
    data_MEG.stim  = data.stim;
    
    %Determine MEG samples corresponding to specific tone
    nTrials  = length(data_MEG.trial);
    nSensors = size( data_MEG.trial{1}, 1 );
    
    samplingFreq    = data_MEG.fsample;
    toneDur_inSecs  = str2num(tonedur_text);
    nSamplesPerTone = toneDur_inSecs * samplingFreq;
    
    nSamplesPerSeries = length(data_MEG.trial{1}(1,:));
    
    t_start_ind(1) = find(data_MEG.time{1} == 0); %Start of first tone defined as t = 0
    t_end_ind(1)   = t_start_ind(1) + nSamplesPerTone -1;
    
    for i_tone = 2:34
        t_start_ind(i_tone) = t_end_ind(i_tone-1) + 1;
        t_end_ind(i_tone)   = t_start_ind(i_tone) + nSamplesPerTone - 1;
    end
    
    %Read out MEG data at specific tone
    %For each trial, extract samples corresponding to respective tone (33), as well as relevant stim/behav info
    cfg = [];
    cfg.begsample = t_start_ind(toneIndex);
    cfg.endsample = t_end_ind(toneIndex);
    data_MEG_Tone33 = ft_redefinetrial(cfg,data_MEG);
    
    %Separate trials according to math. expected tone (p*34)
    filt_low_predp34 = data_MEG.stim.predID == -1;
    filt_middle_predp34 = data_MEG.stim.predID ==  0;
    filt_high_predp34 = data_MEG.stim.predID == 1;
    
    %Select respective trials
    cfg = [];
    cfg.trials = filt_low_predp34;
    data_Tone33_low_predp34 = ft_selectdata(cfg,data_MEG_Tone33);
    
    cfg.trials = filt_middle_predp34;
    data_Tone33_middle_predp34 = ft_selectdata(cfg,data_MEG_Tone33);
    
    cfg.trials = filt_high_predp34;
    data_Tone33_high_predp34 = ft_selectdata(cfg,data_MEG_Tone33);
    
    %Store all trials with identical p*34 in summary struct
    for i_trial = 1:length(data_Tone33_low_predp34.trial)
        dataMEG_p33_allsub.low_predp34{i_sub}(:, :, i_trial) = ...
            data_Tone33_low_predp34.trial{i_trial}(:,:);
    end
    %Output: MEG activity during tone 33 with dims: sensors/timepoints/trials
    
    for i_trial = 1:length(data_Tone33_middle_predp34.trial)
        dataMEG_p33_allsub.middle_predp34{i_sub}(:, :, i_trial) = ...
            data_Tone33_middle_predp34.trial{i_trial}(:,:);
    end
    
    for i_trial = 1:length(data_Tone33_high_predp34.trial)
        dataMEG_p33_allsub.high_predp34{i_sub}(:, :, i_trial) = ...
            data_Tone33_high_predp34.trial{i_trial}(:,:);
    end
    
end

%% 2) Group level averaging (per ToneDur)
%Average across trials
%Output: Matrix with subs*chans*time points
for i_sub = 1:length(subs)
    low_predp34(i_sub,:,:) = ...
        nanmean(dataMEG_p33_allsub.low_predp34{i_sub},3);
    middle_predp34(i_sub,:,:) = ...
        nanmean(dataMEG_p33_allsub.middle_predp34{i_sub},3);
    high_predp34(i_sub,:,:) = ...
        nanmean(dataMEG_p33_allsub.high_predp34{i_sub},3);
end

%Average across subs
%Output: Matrix with Chans*time points
GroupAvg_low_predp34 = squeeze(nanmean(low_predp34,1));
GroupAvg_middle_predp34 = squeeze(nanmean(middle_predp34,1));
GroupAvg_high_predp34 = squeeze(nanmean(high_predp34,1));

SEM_low_predp34 = squeeze(nanstd(low_predp34,1)/sqrt(length(subs))); 
SEM_middle_predp34 = squeeze(nanstd(middle_predp34,1)/sqrt(length(subs)));
SEM_high_predp34 = squeeze(nanstd(high_predp34,1)/sqrt(length(subs)));

%% 3) Define relevant time windows
samplingFreq    = 600;
toneDur_inSecs  = str2num(tonedur_text);
nSamplesPerTone = toneDur_inSecs * samplingFreq;

% define windows
win_size    = 30;
win_overlap = 0;

windows = [1 win_size];
while windows(end,end) < nSamplesPerTone
    windows = [windows; windows(end,:) + (win_size - win_overlap)];
end

if windows(end,end) > nSamplesPerTone
    windows(end,:) = [];
end

windows_inms = (windows / samplingFreq) * 1000;

%% 4) Load in sensor indices for predictive processing clusters (300 ms TD)
%Load in file with predictive processing cluster sensor indices
load([paths_NASTD_MEG.Current_outputdata 'PredEffect/' tonedur_title ...
    '/Group/Group_PredEffect_' tonedur_title '_ClusterCorrect_p0.05_nReps1000.mat']);

%Read out sensor indices per cluster
index_clusterSOI = [];
for i_win = 1:size(windows,1)
    ind_start = windows(i_win, 1);
    ind_end   = windows(i_win, 2);
    
    i_clustercounter = 0;
    for i_cluster = 1:clusters_Exp{i_win}.nClusters
        if clusters_Exp{i_win}.cluster_pval(i_cluster) < 0.025
            index_clusterSOI_proxy= clusters_Exp{i_win}.cluster_sensors{i_cluster};
            
            i_clustercounter = i_clustercounter+1;
            ind_clusterselected_channels{i_win}{i_clustercounter} = index_clusterSOI_proxy;
            index_clusterSOI_proxy = [];
        end
    end
end


%% 5) Place data in plotting-conform structure
for i_win = 1:length(ind_clusterselected_channels)
    for i_cluster = 1:length(ind_clusterselected_channels{i_win})
        if ~isempty(ind_clusterselected_channels{i_win}{i_cluster})
            MEG_GroupAvg_AvgClustChan_lowpredp34trials{i_win}{i_cluster} = ...
                mean(GroupAvg_low_predp34(ind_clusterselected_channels{i_win}{i_cluster},:));
            MEG_GroupAvg_AvgClustChan_middlepredp34trials{i_win}{i_cluster} = ...
                mean(GroupAvg_middle_predp34(ind_clusterselected_channels{i_win}{i_cluster},:));
            MEG_GroupAvg_AvgClustChan_highpredp34trials{i_win}{i_cluster} = ...
                mean(GroupAvg_high_predp34(ind_clusterselected_channels{i_win}{i_cluster},:));
            
            MEG_GroupAvg_SEMClustChan_lowpredp34trials{i_win}{i_cluster} = ...
                mean(SEM_low_predp34(ind_clusterselected_channels{i_win}{i_cluster},:));
            MEG_GroupAvg_SEMClustChan_middlepredp34trials{i_win}{i_cluster} = ...
                mean(SEM_middle_predp34(ind_clusterselected_channels{i_win}{i_cluster},:));
            MEG_GroupAvg_SEMClustChan_highpredp34trials{i_win}{i_cluster} = ...
                mean(SEM_high_predp34(ind_clusterselected_channels{i_win}{i_cluster},:));
        end
    end    
end

%% 6) Plot ERFs per p*34 for selected sensors
if plot_figs == 1
    for i_win = 1:length(ind_clusterselected_channels)
        for i_cluster = 1:length(MEG_GroupAvg_AvgClustChan_lowpredp34trials{i_win})
            if ~isempty(MEG_GroupAvg_AvgClustChan_lowpredp34trials{i_win}{i_cluster})
                
                h = figure;
                set(gcf,'Renderer','painters');
                set(gcf,'units','normalized','outerposition',[0 0 1 1])
                hold on
                
                l1 = shadedErrorBar...
                    (1:length(MEG_GroupAvg_AvgClustChan_lowpredp34trials{i_win}{i_cluster}),...
                    MEG_GroupAvg_AvgClustChan_lowpredp34trials{i_win}{i_cluster}, ...
                    MEG_GroupAvg_SEMClustChan_lowpredp34trials{i_win}{i_cluster});
                    l1.mainLine.Color = [0, 0.4470, 0.7410, 0.7];
                    l1.patch.FaceColor = [0, 0.4470, 0.7410];
                    l1.patch.FaceAlpha = 0.1;
                    l1.mainLine.LineWidth = 4;
                l2 = shadedErrorBar...
                    (1:length(MEG_GroupAvg_AvgClustChan_middlepredp34trials{i_win}{i_cluster}),...
                    MEG_GroupAvg_AvgClustChan_middlepredp34trials{i_win}{i_cluster}, ...
                    MEG_GroupAvg_SEMClustChan_middlepredp34trials{i_win}{i_cluster});
                    l2.mainLine.Color = [0, 0.75, 0.75, 0.7];
                    l2.patch.FaceColor = [0, 0.75, 0.75];
                    l2.patch.FaceAlpha = 0.1;
                    l2.mainLine.LineWidth = 4;
                l3 = shadedErrorBar...
                    (1:length(MEG_GroupAvg_AvgClustChan_highpredp34trials{i_win}{i_cluster}),...
                    MEG_GroupAvg_AvgClustChan_highpredp34trials{i_win}{i_cluster}, ...
                    MEG_GroupAvg_SEMClustChan_highpredp34trials{i_win}{i_cluster});
                    l3.mainLine.Color = [0.8500, 0.3250, 0.0980, 0.7];
                    l3.patch.FaceColor = [0.8500, 0.3250, 0.0980];                    
                    l3.patch.FaceAlpha = 0.1;
                    l3.mainLine.LineWidth = 4;

                    h.CurrentAxes.XLabel.String = 'Time [s]';
                    h.CurrentAxes.XTick = 0:30:length(MEG_GroupAvg_AvgClustChan_lowpredp34trials{i_win}{i_cluster});
                    h.CurrentAxes.XTickLabel = 0:0.05:length(MEG_GroupAvg_AvgClustChan_lowpredp34trials{i_win}{i_cluster})/600;
                    h.CurrentAxes.YLabel.String = 'T';
                
                    xlim([1, length(MEG_GroupAvg_AvgClustChan_lowpredp34trials{i_win}{i_cluster})])
                %             ylim([-0.4e-13 0.6e-13])
                %             ylim([min(all_val)-0.5*max(all_val) max(all_val)+0.5*max(all_val)])
                ylim('auto')
                set(gca,'FontSize',12)
                legend('low p*34', 'middle p*34', 'high p*34')
                
                w1 = num2str( 1000 * windows(i_win, 1) / samplingFreq , 3 );
                w2 = num2str( 1000 * windows(i_win, 2) / samplingFreq , 3 );
                title({['Group-level ERF data during tone 33 for prediction clusterSOI (TW = [' w1 '- ' w2 ' ms]; ClusterNum = ' ...
                    num2str(i_cluster) ', ClusterChan = ' num2str(clusters_Exp{i_win}.cluster_size(i_cluster)) ') as function of p*34']})
                
                if save_figs == 1
                    path_fig = ([paths_NASTD_MEG.Current_outputfig 'ERFp33/' tonedur_title '/Group/']);
                    mkdir(path_fig);
                    figfile = ['Group_ERFp33_perpredp34_Cluster' num2str(i_cluster) 'SOI_TD' tonedur_title '0_TW' w1 '-' w2 'ms.png'];
                    saveas(gcf, [savedir figfile], 'png');
                    close
                end
                
            end            
        end
    end
end

end