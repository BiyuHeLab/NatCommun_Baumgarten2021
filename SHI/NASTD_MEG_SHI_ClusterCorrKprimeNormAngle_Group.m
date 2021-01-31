function NASTD_MEG_SHI_ClusterCorrKprimeNormAngle_Group...
    (param, ...
    paths_NASTD_MEG)
%Aim: Plot ouput of statistical exp vs. shuff comparison for all sensors in
%topoplot - topoplot shows norm/angle and highlights sign. sensors
%(uncorrected p-values for each sensor)

%% 0) Specify vars, paths, and setup fieldtrip
addpath(genpath(paths_NASTD_MEG.ScriptsDir));

%Data input
path_inputdata = [paths_NASTD_MEG.Current_outputdata 'ExpK/TDcomparison/'];
%Data output
path_outputdata = [paths_NASTD_MEG.Current_outputdata 'ExpK/TDcomparison/ClusterData/'];
mkdir(path_outputdata)

%% 1) Load in data
%Data includes Kprime coordinates, vector norm and angle data, and exp vs.
%shuff p-values per sensor and predetermiend sensor-clusters
load([path_inputdata 'Group_KprimeNormAngle.mat']);
load([paths_NASTD_MEG.ScriptsDir 'MEG_sensor_setup_272/label272.mat']); %file with CTF sensor labels for 272 sensors, called 'label'

%% 2) Determine analysis parameters
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

%% 3) Perform cluster-correction for respective input parameter
for i_inputparam = 1:length(param.KprimeComparison.AnalysisParam)
    curr_plotparam = param.KprimeComparison.AnalysisParam{i_inputparam};
    
    for i_win = 1:length(Kprime4AllTD.(curr_plotparam).perSens.exp)
        %Place data (exp input data and uncorrected p values) in topo structure
        topo_stat_exp{i_win} = Kprime4AllTD.(curr_plotparam).perSens.exp{i_win};
        topo_p_exp{i_win}    = Kprime4AllTD.(curr_plotparam).perSens.pval_expvsshuff{i_win};
        
        clusters_exp{i_win} = find_clusters(topo_stat_exp{i_win}, topo_p_exp{i_win}, param.cluster.pval_clusterdef);
        
        %Compute a pseudo-p-value for the shuffled condition used to read out
        %sign. clusters on shuffled data and determine the number of sensors where
        %shuffled values of that repetition are larger than shuffled values from all permutations
        if strcmp(curr_plotparam,'Norm') %Vector norm (p-val determined by shuff > exp)
            for i_sensor = 1:size(Kprime4AllTD.(curr_plotparam).perSens.shuff{i_win},1)
                for i_rep = 1:size(Kprime4AllTD.(curr_plotparam).perSens.shuff{i_win},2)
                    Kprime4AllTD.(curr_plotparam).perSens.pval_shuffvsshuff{i_win}(i_sensor, i_rep) = ...
                        sum(Kprime4AllTD.(curr_plotparam).perSens.shuff{i_win}(i_sensor,:)...
                        >= Kprime4AllTD.(curr_plotparam).perSens.shuff{i_win}(i_sensor,i_rep))...
                        / size(Kprime4AllTD.(curr_plotparam).perSens.shuff{i_win},2);
                end
            end
        else %Vector angle (p-val determined by shuff < exp)
            for i_sensor = 1:size(Kprime4AllTD.(curr_plotparam).perSens.shuff{i_win},1)
                for i_rep = 1:size(Kprime4AllTD.(curr_plotparam).perSens.shuff{i_win},2)
                    Kprime4AllTD.(curr_plotparam).perSens.pval_shuffvsshuff{i_win}(i_sensor, i_rep) = ...
                        sum(Kprime4AllTD.(curr_plotparam).perSens.shuff{i_win}(i_sensor,:)...
                        <= Kprime4AllTD.(curr_plotparam).perSens.shuff{i_win}(i_sensor,i_rep))...
                        / size(Kprime4AllTD.(curr_plotparam).perSens.shuff{i_win},2);
                end
            end            
        end
        
        for i_rep = 1:size(Kprime4AllTD.(curr_plotparam).perSens.shuff{i_win},2)
            topo_stat_shuff{i_rep} = Kprime4AllTD.(curr_plotparam).perSens.shuff{i_win}(:,i_rep);
            topo_p_shuff{i_rep}    = Kprime4AllTD.(curr_plotparam).perSens.pval_shuffvsshuff{i_win}(:,i_rep);
            %find clusters in shuff data
            clusters_shuff{i_rep} = find_clusters(topo_stat_shuff{i_rep}, topo_p_shuff{i_rep}, param.cluster.pval_clusterdef);
            shuffMaxStat{i_win}(i_rep)   = clusters_shuff{i_rep}.maxStatSumAbs;
        end
        
        %Determine sign. sensor-clusters in exp data by testing shuff vs. exp maxstat
        for i_cluster = 1:clusters_exp{i_win}.nClusters
            pval = sum(shuffMaxStat{i_win} >= abs( clusters_exp{i_win}.cluster_statSum(i_cluster) )) / ...
                size(Kprime4AllTD.(curr_plotparam).perSens.shuff{i_win},2);
            clusters_exp{i_win}.cluster_pval(i_cluster) = pval;
        end
        
        %Compute effect size for each cluster
        for i_cluster = 1:clusters_exp{i_win}.nClusters
            filt_cluster = shuffMaxStat{i_win} > 0;
            cluster_effectsize{i_win}{i_cluster} = ...
                (abs(clusters_exp{i_win}.cluster_statSum(i_cluster)) - mean(shuffMaxStat{i_win}(filt_cluster)))...
                / std(shuffMaxStat{i_win}(filt_cluster));
        end
    end
    
    %Save cluster information
    savefile = [path_outputdata 'Group_KprimeClusterCorr_' curr_plotparam '.mat'];
    save(savefile, ...
    'clusters_exp','shuffMaxStat','cluster_effectsize', ... 
    '-v7.3');

%     %Load saved cluster data
%     loadfile = [path_outputdata 'Group_KprimeClusterCorr_' curr_plotparam '.mat'];
%     load(loadfile);
    
    %% 4) Plot resulting clusters as topoplot for respective input parameter    
    %Dtermine z-scaling
    absmax_plot = -Inf;
    absmin_plot = +Inf;
    for i_win = 1:length(Kprime4AllTD.(curr_plotparam).perSens.exp)
        GAvgCluster_topoplot = clusters_exp{i_win}.inputs.topo_stat;
        
        if max(abs(GAvgCluster_topoplot)) > absmax_plot
            absmax_plot = round(max(abs(GAvgCluster_topoplot)),2);
        end
        if min(abs(GAvgCluster_topoplot)) < absmin_plot
            absmin_plot = round(min(abs(GAvgCluster_topoplot)),2);
        end
    end
    
    for i_win = 1:length(Kprime4AllTD.(curr_plotparam).perSens.exp)
        %Determine TW parameters for text plotting
        w1 = num2str(1000 * windows(i_win, 1) / samplingFreq , 3 );
        w2 = num2str(1000 * windows(i_win, 2) / samplingFreq , 3 );
        win_text = [w1 '-' w2 'ms'];
        %Determine if cluster is significant and thus plotted
        n_cluster2plot = 0;
        cluster_sensorindices = [];
        for i_cluster = 1:clusters_exp{i_win}.nClusters%
            cluster_pval = clusters_exp{i_win}.cluster_pval(i_cluster);
            if cluster_pval < param.plot.pval
                n_cluster2plot = n_cluster2plot + 1;
                cluster_sensorindices = [cluster_sensorindices; clusters_exp{i_win}.cluster_sensors{i_cluster}];
                cluster_size(n_cluster2plot) = clusters_exp{i_win}.cluster_size(i_cluster);
            end
        end
        
        %Set up topoplot struct
        dat.dimord = 'chan_time';
        dat.label  = label;
        dat.time   = 0;
        
        cfg = [];
        cfg.layout    = 'CTF275.lay';
        cfg.comment   = 'no';
        cfg.colorbar  = 'yes';
        cfg.zlim      = [absmin_plot absmax_plot];
        cfg.highlight           = 'on';
        cfg.highlightcolor      = [0.9 0.9 0.9];
        cfg.highlightsymbol     = '.';
        cfg.highlightsize       = 30;
        cfg.highlightchannel = dat.label(cluster_sensorindices);
        
        dat.avg = clusters_exp{i_win}.inputs.topo_stat; %perpendicular distance to QLine per sensor
        
        %Plot topoplot
        h = figure;
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        ft_topoplotER(cfg, dat);
        
        %Add title        
        if n_cluster2plot == 0
            cluster_size_text = 'No significant clusters';            
        else
            cluster_size_text = [];
            for i_cluster = 1:n_cluster2plot
                cluster_size_text = [cluster_size_text num2str(cluster_size(i_cluster))];
                if i_cluster < n_cluster2plot
                    cluster_size_text = [cluster_size_text ', '];
                end
            end
        end
        title({['Vector ' curr_plotparam ' - sign. Sensors (exp vs. shuff,p < ' ...
            num2str(param.plot.pval) ' (Cluster-corrected; ' ...
            cluster_size_text ' sign. Sensors); TW: ' win_text]})

        %Save figure
        if param.plot.save   == 1
            path_figs = [paths_NASTD_MEG.Current_outputfig 'ExpKprime/TDcomparison/Topoplot/ClusterCorrected/'];
            mkdir(path_figs)
            
            filename     = ['Group_TopoClusterCorrKprime_' curr_plotparam '_' win_text '.png'];
            figfile      = [path_figs filename];
            saveas(gcf, [figfile], 'png'); %save png version
            delete(h);
        end        
    end
end

end