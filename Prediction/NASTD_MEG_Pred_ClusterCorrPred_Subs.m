function NASTD_MEG_Pred_ClusterCorrPred_Subs ...
    (sub, tonedur_text, ...
    pval_cluster, nReps_cluster, ...
    plot_figs, pval_plotting, save_figs, ...
    paths_NASTD_MEG)
%Aim: Perform cluster correction for single-subject prediction effects.

%% 0) Specify vars, paths, and setup fieldtrip
addpath(genpath(paths_NASTD_MEG.ScriptsDir));
NASTD_MEG_SubInfo

%Tone duration condition for data load-in
if strcmp(tonedur_text,'0.15')
    tonedur_title = '0.15sTD';
elseif  strcmp(tonedur_text,'0.3')
    tonedur_title = '0.3sTD';
elseif strcmp(tonedur_text,'0.6')
    tonedur_title = '0.6sTD';
end

%Data input
path_input = [paths_NASTD_MEG.Current_outputdata 'PredEffect/' tonedur_title '/' sub '/'];
%Data output
path_output = [paths_NASTD_MEG.Current_outputdata 'PredEffect/' tonedur_title '/' sub '/'];

%% 1) Load and prepare input data
load([path_input sub '_PredEffect_' tonedur_title '.mat']);

%% 2) Define analysis parameters
win_size    = 30;
win_overlap = 0;

samplingFreq = 600;
toneDur_inSecs  = str2num(tonedur_text);
nSamplesPerTone = toneDur_inSecs * samplingFreq;

windows = [1 win_size];
while windows(end,end) < nSamplesPerTone
    windows = [windows; windows(end,:) + (win_size - win_overlap)];
end

if windows(end,end) > nSamplesPerTone
    windows(end,:) = [];
end

windows_inms = (windows / samplingFreq) * 1000;
nWindows = size(windows,1);

%% 3) Prepare input for cluster analysis
%Read out dependent (MEG activity at tone 33) and independent variable (p*34)
dv_temp = ERF_p33_win; %dependent var: MEG activity at tone 33
iv_temp = series_predp34_discretized; %independent var: p*34

%Change format to timewin, sensor, trial (i.e., seperate sensors into diff cells)
nSensors = size(dv_temp{1}, 1);
for i_win = 1:nWindows
    for i_sensor = 1:nSensors
        dv_forCluster_win{i_win}{i_sensor} = dv_temp{i_win}(i_sensor, :)';
    end
end

iv_forCluster = iv_temp';

clear dv_temp iv_temp

%% 4) Perform cluster correction
%Aim: Perform linear regression between DV and IV i for each sensor.
%The resulting t-statistic is used as parameter entering cluster analysis.
%Clusters are computed for Exp and Shuffled distibutions.
%For shuffled, the DV trialorder is shuffled and for each permutation
%paired with the original IV in order to create a null distribution.

%Output; Stat parameter (per sensor) for exp distribution is the t-stat from the
%t-test compairing Beta value of linear regression across subjects. p values
%are determined by as the ratio of Exp max cluster stat > Shuffled max cluster stat.
for i_win = 1:nWindows
    dv_forCluster = dv_forCluster_win{i_win};
    [clusters_Exp{i_win}, clusters_Shuff{i_win}, clusters_Shuff_MaxStat{i_win}] = ...
        NASTD_MEG_Pred_ClusterCorrect_Subs(dv_forCluster, iv_forCluster, pval_cluster, nReps_cluster);
end

%% 5) Save output
outputfile_name = [sub '_PredEffect_' tonedur_title ...
    '_ClusterCorrect_p' num2str(pval_cluster) ...
    '_nReps' num2str(nReps_cluster) '.mat'];

save([path_output outputfile_name], ...
    'clusters_Exp', 'clusters_Shuff_MaxStat',...
    '-v7.3');

%% 6) Plot topoplot (showing cluster-corrected p-values per sensor for each time window)
if plot_figs == 1
    %Load in label file
    load([paths_NASTD_MEG.ScriptsDir 'MEG_sensor_setup_272/label272.mat'] ); %file with CTF sensor labels for 272 sensors
    
    %Determine z-scaling (stable across time windows)
    %Check if there are any sign cluster present for time window
    absmax_plot = 0;
    for i_win = 1:nWindows
        nCluster2plot = 0;
        
        for i_cluster = 1:clusters_Exp{i_win}.nClusters %loop across clusters present in time window
            
            cluster_p = clusters_Exp{i_win}.cluster_pval(i_cluster); %read out pval for current cluster
            if cluster_p < pval_plotting %Check if pval is below plotting threshold
                nCluster2plot = nCluster2plot + 1; %If yes, add cluster to plot list
                
                Cluster_SensIndices  = clusters_Exp{i_win}.cluster_sensors{i_cluster}; %Collect all sign. sensors in one file
                
                topo_plot(Cluster_SensIndices) = clusters_Exp{i_win}.inputs.topo_stat(Cluster_SensIndices); %Read out stats for said sensors
                Cluster_NumSens(nCluster2plot) = clusters_Exp{i_win}.cluster_size(i_cluster); %Read out number of sensors for each cluster
                
            end
        end
        
        %Read out maximum stat size to determine upper limit for scaling
        if nCluster2plot > 0
            if max(abs(topo_plot)) > absmax_plot
                absmax_plot = max(abs(topo_plot));
            end
        end
        
        %Update this for each window to find upper limit for scaling across windows
        topo_plot = clusters_Exp{i_win}.inputs.topo_stat;
        if max(abs(topo_plot)) > absmax_plot
            absmax_plot = max(abs(topo_plot));
        end
        
    end
    absmax_plot = round(absmax_plot);
    
    %Create topo struct with group-avg prediction effects per sensor and sign. cluster sensor-groups
    for i_win = 1:nWindows
        
        %Set up title information
        w1 = num2str( 1000 * windows(i_win, 1) / samplingFreq , 3 );
        w2 = num2str( 1000 * windows(i_win, 2) / samplingFreq , 3 );
        win_title = ['time window = [' w1 ' ms - ' w2 ' ms]'];
        
        assoc_title = 'Prediction effect (t-stat from linear regression)';
        title_text = ['evaluated at tone #' num2str(toneIndex) ...
            ', prediction using k = ' num2str(predictive_sequencerange) ' tones'];
                
        %Read out indices and number of sensors per cluster
        nSensors = length(clusters_Exp{1}.topo_cluster);
        topo_plot = clusters_Exp{i_win}.inputs.topo_stat; %Copy prediction stats for all sens
        
        nClusterPlot = 0;
        Cluster_SensIndices = [];
        
        for i_cluster = 1:clusters_Exp{i_win}.nClusters
            
            cluster_p(i_cluster) = clusters_Exp{i_win}.cluster_pval(i_cluster);
            if cluster_p(i_cluster) < pval_plotting %if cluster is sign.
                nClusterPlot = nClusterPlot + 1; %add cluster to plot list
                
                Cluster_SensIndices = [Cluster_SensIndices clusters_Exp{i_win}.cluster_sensors{i_cluster}];
                
                Cluster_NumSens(nClusterPlot) = clusters_Exp{i_win}.cluster_size(i_cluster);
                
            end            
        end
        
        clear dat
        dat.dimord = 'chan_time';
        dat.label  = label;
        dat.time   = 0;
        dat.avg = topo_plot';
        
        %3.3 Plot topo
        cfg = [];
        cfg.layout    = 'CTF275.lay';
        cfg.colorbar  = 'yes';
        cfg.comment   = 'no';
        cfg.zlim = [-absmax_plot, absmax_plot];
        
        cfg.highlight        = 'on';
        cfg.highlightchannel = dat.label(Cluster_SensIndices);
        cfg.highlightsymbol  = '.';
        cfg.highlightcolor   = [1 1 1];
        cfg.highlightsize    = 25;
        
        h=figure;
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        
        ft_topoplotER(cfg, dat);
                
        if nClusterPlot == 0
            cluster_size_text = 'No significant clusters';
        else
            cluster_size_text = [];
            for i_cluster = 1:nClusterPlot
                cluster_size_text = [cluster_size_text num2str(Cluster_NumSens(i_cluster))];
                if i_cluster < nClusterPlot
                    cluster_size_text = [cluster_size_text ', '];
                end
            end
        end
        
        title({[sub ' ' assoc_title ' p < ' num2str(pval_plotting) ...
            ', cluster-corrected; cluster sizes = ' cluster_size_text], ['TD: ' tonedur_title '; ERF ' win_title], title_text})
        
        %Save topoplot
        if save_figs    == 1
            path_fig = ([paths_NASTD_MEG.Current_outputfig 'PredEffect_uncorrected/' tonedur_title '/' sub '/']);
            mkdir(path_fig);
            filename     = [sub '_PredEffect_ClusterCorr_' tonedur_title '_TW' w1 '-' w2 '.png'];
            figfile      = [path_fig filename];
            saveas(gcf, [figfile], 'png'); 
            delete(h);
        end
        
    end
end

end