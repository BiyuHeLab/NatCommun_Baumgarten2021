function NASTD_MEG_ERF_StatM100SpatialFilter...
    (i_tonedur, ...
    param, paths_NASTD_MEG)

%Aim: Compute group-level statistics for comparison between low vs. high p*34 
%for M100 sensory filter weighted data

%% 1) Prepare input (single-subject M100 sensory filter weighted data)
path_input = [paths_NASTD_MEG.Current_outputdata 'Data_M100SpatialFilter/Ssub/' param.tonedur_text{i_tonedur} 'sTD/'];
path_fig = [paths_NASTD_MEG.Current_outputfig 'M100SpatialFilter/GroupStat/' param.tonedur_text{i_tonedur} 'sTD/'];
path_output = [paths_NASTD_MEG.Current_outputdata 'Data_M100SpatialFilter/GroupStat/' param.tonedur_text{i_tonedur} 'sTD/'];
mkdir([path_output]);
mkdir(path_fig)

%Select subcondition trials from individual M100 spatially filtered data
sub_counter = 0;
for i_sub = param.M100.Subs{i_tonedur}
    
    sub_counter = sub_counter+1;  
    
    %Load in weighted single-trial data
    load([path_input  param.subs{i_sub} '_M100SpatFilt' ...
        '_WeightedData_' param.tonedur_text{i_tonedur} 'sTD.mat']);    
    
    %Determine tone-dependent begin/end time and time window for stat comparison
    if i_sub == 1
        nSamplesPerTone = str2double(param.tonedur_text{i_tonedur}) * param.samplefreq;
        index_T1start = find(ERF_M100filtweight_perTrial.time{1} == 0); %determine start of first tone
        Samples_AllToneStart = [index_T1start:nSamplesPerTone:(index_T1start+34*nSamplesPerTone)];
        trial_time = ERF_M100filtweight_perTrial.time{1}; %Save trial time date in seperate file
        
        TOI_statcomp = [ERF_M100filtweight_perTrial.time{1}(Samples_AllToneStart(1)) ...
            ERF_M100filtweight_perTrial.time{1}(Samples_AllToneStart(35))+0.4];
            %Start p1 until start resp window
    end
    
    %Select conditions of interest
    %low p*34
    cfg                 = [];
    cfg.trials          = ERF_M100filtweight_perTrial.trialindices.predp34_low;
    ERF_M100filtweight_Predp34Low{sub_counter}    = ft_preprocessing(cfg, ERF_M100filtweight_perTrial); %Select trials
    cfg                 = [];
    cfg.keeptrials      = 'no';
    ERF_M100filtweight_Predp34Low{sub_counter}    = ft_timelockanalysis(cfg,ERF_M100filtweight_Predp34Low{sub_counter}); %Avg over trials
    cfg = [];
    cfg.avgoverchan = 'yes'; %Avg over channels
    ERF_M100filtweight_Predp34Low{sub_counter} = ft_selectdata(cfg,ERF_M100filtweight_Predp34Low{sub_counter});
        
    %high p*34
    cfg                 = [];
    cfg.trials       	= ERF_M100filtweight_perTrial.trialindices.predp34_high;
    ERF_M100filtweight_Predp34High{sub_counter} 	= ft_preprocessing(cfg, ERF_M100filtweight_perTrial);
    cfg                 = [];
    cfg.keeptrials      = 'no';
    ERF_M100filtweight_Predp34High{sub_counter} = ft_timelockanalysis(cfg,ERF_M100filtweight_Predp34High{sub_counter});
    cfg = [];
    cfg.avgoverchan = 'yes'; %Avg over channels
    ERF_M100filtweight_Predp34High{sub_counter} = ft_selectdata(cfg,ERF_M100filtweight_Predp34High{sub_counter});
    
    %Place data in common struct
    Ssub_SpatFiltData.Predp34Low(:,:,sub_counter)    = ERF_M100filtweight_Predp34Low{sub_counter}.avg;
    Ssub_SpatFiltData.Predp34High(:,:,sub_counter)   = ERF_M100filtweight_Predp34High{sub_counter}.avg;
    
    clear ERF_M100filtweight_perTrial
end

%% 2) Average data across subjects
%Non-FT struct
%Average across channels for each sub
Ssub_SensAvg_SpatFiltData.Predp34Low = squeeze(mean(Ssub_SpatFiltData.Predp34Low,1));
Ssub_SensAvg_SpatFiltData.Predp34High = squeeze(mean(Ssub_SpatFiltData.Predp34High,1));

%Average across subejcts for each sample
GAvg_SensAvg_SpatFiltData.Predp34Low = mean(Ssub_SensAvg_SpatFiltData.Predp34Low,2);
GAvg_SensAvg_SpatFiltData.Predp34High = mean(Ssub_SensAvg_SpatFiltData.Predp34High,2);

%% 3) Perform Group-level cluster-statistics
   
        cfg                     = [];
        cfg.channel             = 'all'; %1 spatial filter weighted channel
        cfg.latency             = TOI_statcomp;
        cfg.avgovertime         = 'no';
        cfg.parameter           = 'avg';

        cfg.method              = 'montecarlo';
        cfg.alpha               = param.stat.pval_clusterstat;
        cfg.tail                = 0;        %two-sided test
        cfg.correcttail         = 'prob';   %correct prob.subfield for two-sided test

        cfg.correctm            = 'cluster';
        cfg.clusteralpha        = 0.05;
        cfg.clusterstatistic    = 'maxsum';
        cfg.clustertail         = 0;

        cfg.numrandomization    = param.stat.NumReps;

        Nsub                    = length(ERF_M100filtweight_Predp34Low);
        cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
        cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

        cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
        cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
        
        cfg.statistic           = 'ft_statfun_depsamplesT'; %within-subjects comparison

        StatGAvg_Predp34LowHigh.stat = ft_timelockstatistics(cfg,...
                ERF_M100filtweight_Predp34Low{:},...
                ERF_M100filtweight_Predp34High{:});
            
%) Save stat-data
    savefile = [path_output  'GroupClustStatData_n' num2str(length(ERF_M100filtweight_Predp34Low)) ...
        '_M100SpatFilt_' param.tonedur_text{i_tonedur} 'TD.mat'];
    save(savefile, ...
        'nSamplesPerTone', 'Samples_AllToneStart', 'trial_time', ...
        'GAvg_SensAvg_SpatFiltData',  'StatGAvg_Predp34LowHigh', ...
        'param', '-v7.3');
           
%% 4) Plot  comparison between conditions (with stats)
% load([path_output  'GroupClustStatData_n' num2str(length(param.M100.Subs{i_tonedur})) ...
%         '_M100SpatFilt_' param.tonedur_text{i_tonedur} 'TD.mat']);

%Set up figure
f3 = figure;
set(gcf,'Renderer','painters');
set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen

%Determine Y axis limits
maxY = max([GAvg_SensAvg_SpatFiltData.Predp34Low; ...
    GAvg_SensAvg_SpatFiltData.Predp34High]);
minY = min([GAvg_SensAvg_SpatFiltData.Predp34Low; ...
    GAvg_SensAvg_SpatFiltData.Predp34High]);

%find sample when FTPL ratings starts - we only want to plot until there
[~, sample_FTPLrating] = min(abs(trial_time - ...
    (trial_time(Samples_AllToneStart(34)+nSamplesPerTone)+0.4)));
TP_FTPLrating = trial_time(sample_FTPLrating);

%Plot grandaverage traces for low vs. high p*34
hold on;
l1 = plot(1:length(GAvg_SensAvg_SpatFiltData.Predp34Low(1:sample_FTPLrating)),...
    GAvg_SensAvg_SpatFiltData.Predp34Low(1:sample_FTPLrating));
l1.Color = [0, 0.4470, 0.7410];
l1.LineWidth = 3;
l3 = plot(1:length(GAvg_SensAvg_SpatFiltData.Predp34High(1:sample_FTPLrating)),...
    GAvg_SensAvg_SpatFiltData.Predp34High(1:sample_FTPLrating));
l3.Color = [0.9290, 0.6940, 0.1250];
l3.LineWidth = 3;

%Highlights samples showing sign. differences
sign_samples = find(StatGAvg_Predp34LowHigh.stat.mask); %highlight sign. samples
%CAVE: Stat has different sample-to-timeppoint allocation - adjust
samplediff = find(StatGAvg_Predp34LowHigh.stat.time(1) == trial_time);
adjusted_sign_samples = zeros(1,length(trial_time));
adjusted_sign_samples(sign_samples + (samplediff)) = 1;

hold on;
area(1:length(trial_time(1:sample_FTPLrating)),...
    [adjusted_sign_samples(1:sample_FTPLrating) ],...
    'basevalue',0,'FaceColor',[0.69, 0, 0.16],...
    'FaceAlpha', 0.5,'LineStyle','none');
hold on;
area(1:length(trial_time(1:sample_FTPLrating)),...
    -adjusted_sign_samples(1:sample_FTPLrating),...
    'basevalue',0,'FaceColor',[0.69, 0, 0.16],...
    'FaceAlpha', 0.5,'LineStyle','none');

%plot vertical lines for tone starts
for i_sample = 1:length(Samples_AllToneStart)
    hold on;
    plot([Samples_AllToneStart(i_sample) Samples_AllToneStart(i_sample)] , ...
        [f3.CurrentAxes.YLim(1) f3.CurrentAxes.YLim(2)],'Color',[0 0 0 0.5],'LineWidth',0.5)
end


%plot horizontal 0 line
hold on;
plot([1 length(GAvg_SensAvg_SpatFiltData.Predp34Low(1:sample_FTPLrating))] , ...
    [0 0],'Color',[0 0 0 0.5],'LineWidth',0.5)

%Add legend and axis labels
legend('Low p*34', 'High p*34','p < 0.05 cluster-corrected')
f3.CurrentAxes.XLim = [1 sample_FTPLrating];
Samples_AllToneStart = [index_T1start : nSamplesPerTone : (index_T1start + (34*nSamplesPerTone))];
Samples_StartRespDisplay = Samples_AllToneStart(end) + (param.samplefreq*0.4); %400ms poststim time
Samples_AllToneStart = [Samples_AllToneStart, Samples_StartRespDisplay];
f3.CurrentAxes.XTick = Samples_AllToneStart;
f3.CurrentAxes.XTickLabel = {'p1','','','','','','','','','p10',...
    '','','','','','','','','','p20',...
    '','','','','','','','','','p30',...
    '','','p33','p34','','RespDisplay'};
f3.CurrentAxes.YLabel.String = 'Amplitude (weighted avg across sensors)';
f3.CurrentAxes.YLim = [minY*1.1 maxY*1.1];

title({['Group-level Stat. Comparison (Cluster-corrected; n = ' num2str(length(param.M100.Subs{i_tonedur})) ') - ' ...
    param.tonedur_text{i_tonedur} ' TD'] ...
    ['Neuromagnetic activity low vs. high p*34'] ...
    ['M100 Spatial Filter (averaged across all (weighted) sensors)']});

if param.plot.save
    Figtitle = (['GroupClustStatData_n' num2str(length(param.M100.Subs{i_tonedur})) '_' param.tonedur_text{i_tonedur} ...
        'sTD_SpatialFilterM100_' num2str(param.M100.TW(1)) 'to' num2str(param.M100.TW(2))...
        's.png']);
    saveas(gcf, [path_fig Figtitle], 'png'); %save png version
    close
end

end