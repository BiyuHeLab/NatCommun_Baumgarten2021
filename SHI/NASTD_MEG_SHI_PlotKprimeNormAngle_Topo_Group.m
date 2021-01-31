function NASTD_MEG_SHI_PlotKprimeNormAngle_Topo_Group...
    (param, ...
    paths_NASTD_MEG)
%Aim: Plot ouput of statistical exp vs. shuff comparison for all sensors in
%topoplot - topoplot shows norm/angle and highlights sign. sensors
%(uncorrected p-values for each sensor)

%% 0) Specify vars, paths, and setup fieldtrip
addpath(genpath(paths_NASTD_MEG.ScriptsDir));

%Data output
path_inputdata = [paths_NASTD_MEG.Current_outputdata 'ExpK/TDcomparison/'];

%% 1) Load in data
%Data includes Kprime coordinates, vector norm and angle data, and exp vs.
%shuff p-values per sensor and predetermiend sensor-clusters
load([path_inputdata 'Group_KprimeNormAngle.mat']);
load([paths_NASTD_MEG.ScriptsDir 'MEG_sensor_setup_272/label272.mat']); %file with CTF sensor labels for 272 sensors, called 'label'

%% 2) Plot topoplot for respective parameter
%Define time windows used for analysis
win_size    = 30;
win_overlap = 0;
samplingFreq = 600;
nSamplesPerTone = 0.15 * samplingFreq; %shortest tone dur * samplingFreq
%Define number, start and end sample of window per tone
windows = [1 win_size];
while windows(end,end) < nSamplesPerTone
    windows = [windows; windows(end,:) + (win_size - win_overlap)];
end
if windows(end,end) > nSamplesPerTone
    windows(end,:) = [];
end
%Compute window start/end time in ms for each time window
windows_inms = (windows / samplingFreq) * 1000;

for i_plotparam = 1:length(param.KprimeComparison.AnalysisParam)
    curr_plotparam = param.KprimeComparison.AnalysisParam{i_plotparam};
    
    %Determine color scaling (consistent across TW)
    for i_win = 1:length(Kprime4AllTD.Norm.perSens.exp)
        maxZ(i_win,1) = ceil(max(Kprime4AllTD.(curr_plotparam).perSens.exp{i_win}));
        minZ(i_win,1) = floor(min(Kprime4AllTD.(curr_plotparam).perSens.exp{i_win}));
    end
    maxZ = max(maxZ);
    minZ = min(minZ);
    
    for i_win = 1:length(Kprime4AllTD.(curr_plotparam).perSens.exp)        
        %Determine TW parameters for text plotting
        w1 = num2str(1000 * windows(i_win, 1) / samplingFreq , 3 );
        w2 = num2str(1000 * windows(i_win, 2) / samplingFreq , 3 );
        win_text = [w1 '-' w2 'ms'];
        
        %Set up topoplot struct
        dat.dimord = 'chan_time';
        dat.label  = label;
        dat.time   = 0;
        
        cfg = [];
        cfg.layout    = 'CTF275.lay';
        cfg.comment   = 'no';
        cfg.colorbar  = 'yes';
        cfg.zlim      = [minZ maxZ];
        cfg.highlight           = 'on';
        cfg.highlightcolor      = [0.9 0.9 0.9];
        cfg.highlightsymbol     = '.';
        cfg.highlightsize       = 30;
        cfg.highlightchannel = find(Kprime4AllTD.(curr_plotparam).perSens.pval_expvsshuff{i_win} < param.plot.pval);
        
        dat.avg = Kprime4AllTD.(curr_plotparam).perSens.exp{i_win}; %perpendicular distance to QLine per sensor
        
        %Plot topoplot
        h = figure;
        set(gcf,'units','normalized','outerposition',[0 0 1 1])        
        ft_topoplotER(cfg, dat);
        
        %Add title
        title({['Vector ' curr_plotparam ' - sign. Sensors (exp vs. shuff,p < ' num2str(param.plot.pval) ' (uncorrected); TW: ' win_text]})
        
        %Save figure
        if param.plot.save   == 1
            path_figs = [paths_NASTD_MEG.Current_outputfig 'ExpKprime/TDcomparison/Topoplot/Uncorrected/'];
            mkdir(path_figs)
            
            filename     = ['Group_TopoUncorrectedKprime_' curr_plotparam '_' win_text '.png'];
            figfile      = [path_figs filename];
            saveas(gcf, [figfile], 'png'); %save png version
            delete(h);
        end
    end
end

end