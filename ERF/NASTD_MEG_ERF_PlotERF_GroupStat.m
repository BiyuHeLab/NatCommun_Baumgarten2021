function NASTD_MEG_ERF_PlotERF_GroupStat...
    (ERF_GroupStat_LowvsHighPredp34, i_tonedur, param, paths_NASTD_MEG)

%Aim: Plot ERF across the whole tone sequence for low vs. high p*34 and
%highlight samples showing sign. difference between conditions

%Set up figure
    figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1])

%Plot group-level traces of both conditions
    plot(1:length(ERF_GroupStat_LowvsHighPredp34.avg.predp34_low.avg), ...
        ERF_GroupStat_LowvsHighPredp34.avg.predp34_low.avg,'b','LineWidth',3);
    hold on; 
    plot(1:length(ERF_GroupStat_LowvsHighPredp34.avg.predp34_high.avg), ...
        ERF_GroupStat_LowvsHighPredp34.avg.predp34_high.avg,'Color',[0.9290, 0.6940, 0.1250],'LineWidth',3);    
        
    %Determine samples corresponding with tone starts and plot them on x-axis
    StartSample_Tone = ...
        zeros(1,length(ERF_GroupStat_LowvsHighPredp34.avg.predp34_low.avg));
    nSamplesPerTone = str2double(param.tonedur_text{i_tonedur}) * param.samplefreq;
    index_T1start = find(round(ERF_GroupStat_LowvsHighPredp34.avg.predp34_low.time,4) == 0);    
    Samples_AllToneStart = [index_T1start : nSamplesPerTone : (index_T1start + (34*nSamplesPerTone))];
    Samples_StartRespDisplay = Samples_AllToneStart(end) + (param.samplefreq*0.4); %400ms poststim time
    Samples_AllToneStart = [Samples_AllToneStart, Samples_StartRespDisplay];
    Label_AllToneStart = {'p1','','','','','','','','','p10',...
        '','','','','','','','','','p20',...
        '','','','','','','','','','p30',...
        '','','p33','p34','','RespDisplay'};
    
    for i_SampleToneStart = Samples_AllToneStart
        StartSample_Tone(i_SampleToneStart) = 1;
    end
    
    plot(1:length(StartSample_Tone),...
        StartSample_Tone,...
        'k','LineWidth',0.1) %plot mark for each tone-start-sample
    plot(1:length(StartSample_Tone),...
        StartSample_Tone * (-1),...
        'k','LineWidth',0.1) %plot mark for each tone-start-sample
    
    %Shade significant samples
    sign_samples = find(ERF_GroupStat_LowvsHighPredp34.stat_LowvsHighPredp34.mask); %read out sign. samples
    %CAVE: .stat subfield has different sample-to-timeppoint allocation because
    %statistical analysis did not take prestimulus duration into account - adjust
    samplediff = ...
        find(ERF_GroupStat_LowvsHighPredp34.stat_LowvsHighPredp34.time(1) == ...
        ERF_GroupStat_LowvsHighPredp34.avg.predp34_low.time);
    adjusted_sign_samples = zeros(1,length(ERF_GroupStat_LowvsHighPredp34.avg.predp34_low.time));
    adjusted_sign_samples(sign_samples + (samplediff)) = 1;
        
    grey  = [127 127 127]./255;
    area(1:length(adjusted_sign_samples), ...
        adjusted_sign_samples,'basevalue',0,'FaceColor',grey,'FaceAlpha', 0.3, 'LineStyle','none');
    area(1:length(adjusted_sign_samples), ...
        adjusted_sign_samples * (-1),'basevalue',0,'FaceColor',grey,'FaceAlpha', 0.3, 'LineStyle','none');

    %Determine axes
    xlim([0 length(ERF_GroupStat_LowvsHighPredp34.avg.predp34_low.avg ...
        (1:Samples_StartRespDisplay))]) %%-0.5 to start resp window
    ylim([min([min(ERF_GroupStat_LowvsHighPredp34.avg.predp34_low.avg) ... %max min)
        min(ERF_GroupStat_LowvsHighPredp34.avg.predp34_high.avg)]) ...
        max([max(ERF_GroupStat_LowvsHighPredp34.avg.predp34_low.avg) ...
        max(ERF_GroupStat_LowvsHighPredp34.avg.predp34_high.avg)])])    
    
    %x-axis shows tones [tone number]
    set(gca,'FontSize',10,'XTickLabelRotation',45)
    set(gca,'xtick',Samples_AllToneStart)
    x_label = Label_AllToneStart;
    set(gca,'FontSize',10)
    xticklabels(x_label);
    
    %Add legend
    conditions = fieldnames(ERF_GroupStat_LowvsHighPredp34.avg);
    l = legend(conditions{1}, conditions{2}, 'Location', 'best');
    set(l,'Interpreter', 'none');
    
    %Add title
    title({['Group-level ERF incl. stat. comparison (cluster-corrected) for low vs. high p*34;  - Tone Dur: ' param.tonedur_text{i_tonedur}]});
    
    %Save fig
    if param.plot.save == 1        
        path_fig = ([paths_NASTD_MEG.Current_outputfig 'GAvg/Stat/LowvsHighPredp34/']);
        mkdir(path_fig);
        filename = ['GroupLevelERF_StatLowvsHighPredp34_' ...
            param.tonedur_text{i_tonedur} 'sTD.png'];
        figfile      = [path_fig filename];
        saveas(gcf, [figfile], 'png'); %save png version
        close
    end
    
end