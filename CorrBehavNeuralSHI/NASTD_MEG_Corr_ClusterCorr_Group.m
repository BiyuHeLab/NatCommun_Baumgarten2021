function NASTD_MEG_Corr_ClusterCorr_Group ...
    (Kprime_AllSub, Fstatinteraction_FTPLrating, ...
    rho_CorrexpKFstat_allSens, pval_CorrexpKFstat_allSens, ...
    param_NASTD_MEG, ...
    paths_NASTD_MEG)

%Aim: Correlate Kprime and Fstatistics for each sensor, read out and plot
%significant cluster-corrected sensor-clusters.

%% 1) Determine parameters
%Determine Window length
win_size    = 30;
win_overlap = 0;
samplingFreq = 600;
toneDur_inSecs = 0.15;
nSamplesPerTone = toneDur_inSecs * samplingFreq;
windows = [1 win_size];

while windows(end,end) < nSamplesPerTone
    windows = [windows; windows(end,:) + (win_size - win_overlap)];
end
if windows(end,end) > nSamplesPerTone
    windows(end,:) = [];
end
windows_inms = (windows / samplingFreq) * 1000;

%% 2. Construct null distribution by correlating behavioral performance with shuffled K-values
%Correlate behavioral performance with a random selection of 1 of the 100
%shuffled runs for each subject. Repeat this 1000 times (with replacement)
%to receive null distribution.

%Create shuffled null distribution based on shuffled Kprime values (averaged across TD)
%by averaging across tone durs within each rep within each subject
for i_win = 1:length(Kprime_AllSub.Shuff{1})
    for i_sub = 1:size(Kprime_AllSub.Shuff,1)
        for i_rep = 1:size(Kprime_AllSub.Shuff{i_sub,1}{i_win},2)
            
            Kprime_AllSub.Shuff{i_sub,4}{i_win}(:,i_rep) = ...
                mean([Kprime_AllSub.Shuff{i_sub,1}{i_win}(:,i_rep), ...
                Kprime_AllSub.Shuff{i_sub,2}{i_win}(:,i_rep), ...
                Kprime_AllSub.Shuff{i_sub,3}{i_win}(:,i_rep)],2);
            
        end
    end
end

%Create shuffled K-prime null distribution (random selection of 1
%shuffled K-prime (K-prime per sensor) run per subject, repated 1000 times
%(with replacement)
for i_win = 1:length(Kprime_AllSub.Shuff{1,4})
    for i_sub = 1:size(Kprime_AllSub.Shuff,1)
        for i_rep = 1:param_NASTD_MEG.cluster.NumReps
            Kprime_ShuffNullDist{i_win}{i_rep}(i_sub,:) = ...
                Kprime_AllSub.Shuff{i_sub,4}{i_win} ...
                (:,randi([1 size(Kprime_AllSub.Shuff{1,4}{i_win},2)]));
        end
    end
    
    %Correlate Kprime and F-statistics across subjects (Per sensor,for all sensors)
    for i_rep = 1:param_NASTD_MEG.cluster.NumReps
        [rho_CorrshuffKFstat_persens{i_win}(:,i_rep), ...
            pval_CorrshuffKFstat_persens{i_win}(:,i_rep)] = ...
            corr(Kprime_ShuffNullDist{i_win}{i_rep}, ...
            Fstatinteraction_FTPLrating(:), ...
            'type', 'Spearman');
    end
    
    %% 3. Perform cluster analysis
    %Place data (rho and p values) in topo form and read out sensor-clusters for exp data
    topo_stat_exp{i_win} = rho_CorrexpKFstat_allSens{i_win};
    topo_p_exp{i_win}    = pval_CorrexpKFstat_allSens{i_win};
    
    cluster_exp{i_win} = find_clusters...
        (topo_stat_exp{i_win}, topo_p_exp{i_win}, param_NASTD_MEG.cluster.pval_clusterdef);
    
    %Place data (rho and p values) in topo form and read out sensor-clusters for shuff data
    for i_rep = 1:param_NASTD_MEG.cluster.NumReps
        topo_stat_shuff{i_rep} = rho_CorrshuffKFstat_persens{i_win}(:,i_rep);
        topo_p_shuff{i_rep}    = pval_CorrshuffKFstat_persens{i_win}(:,i_rep);
        
        %find clusters in shuff data
        clusters_shuff{i_rep} = find_clusters...
            (topo_stat_shuff{i_rep}, topo_p_shuff{i_rep}, param_NASTD_MEG.cluster.pval_clusterdef);
        shuffMaxStat{i_win}(i_rep)   = clusters_shuff{i_rep}.maxStatSumAbs;
    end
    
    %Determine sign. sensor-clusters in exp data by testing shuff vs. exp maxstat
    for i_clust = 1:cluster_exp{i_win}.nClusters
        pval = sum(shuffMaxStat{i_win} >= abs( cluster_exp{i_win}.cluster_statSum(i_clust) )) ...
            / param_NASTD_MEG.cluster.NumReps;
        cluster_exp{i_win}.cluster_pval(i_clust) = pval;
    end
    
    %Compute effect size for each cluster
    for i_cluster = 1:cluster_exp{i_win}.nClusters
        filt_cluster = shuffMaxStat{i_win} > 0;
        cluster_effectsize{i_win}{i_cluster} = ...
            (abs(cluster_exp{i_win}.cluster_statSum(i_cluster)) - ...
            mean(shuffMaxStat{i_win}(filt_cluster)))...
            / std(shuffMaxStat{i_win}(filt_cluster));
    end
end

%% 4. Save cluster data
%Save cluster information
% path_outputdata = [paths_NASTD_MEG.Current_outputdata 'ClusterCorr/'];
% mkdir(path_outputdata)
% savefile = [path_outputdata 'Group_CorrKprimeRho_ClusterData.mat'];
% save(savefile, ...
%     'cluster_exp','shuffMaxStat','cluster_effectsize', ...
%     '-v7.3');

%% 5. Plot significant clusters in topoplot
if param_NASTD_MEG.plot.plot == 1
    load([paths_NASTD_MEG.ScriptsDir 'MEG_sensor_setup_272/label272.mat']); %file with CTF sensor labels for 272 sensors, called 'label'
    
    for i_win = 1:length(cluster_exp)
        
        %Create topo struct
        w1 = num2str( 1000 * windows(i_win, 1) / samplingFreq , 3 );
        w2 = num2str( 1000 * windows(i_win, 2) / samplingFreq , 3 );
        
        nSensors = length(cluster_exp{i_win}.topo_cluster);
        data_topoplot = cluster_exp{i_win}.inputs.topo_stat; %Group-avg Exp k-prime vals per sensor
        
        clear dat
        dat.dimord = 'chan_time';
        dat.label  = label;
        dat.time   = 0;
        dat.avg = data_topoplot;
        
        %Determine if cluster is significant and thus plotted
        n_cluster2plot = 0;
        cluster_sensorindices = [];
        for i_cluster = 1:cluster_exp{i_win}.nClusters
            cluster_pval = cluster_exp{i_win}.cluster_pval(i_cluster);
            if cluster_pval < param_NASTD_MEG.plot.pval
                n_cluster2plot = n_cluster2plot + 1;
                cluster_sensorindices = [cluster_sensorindices; cluster_exp{i_win}.cluster_sensors{i_cluster}];
                cluster_size(n_cluster2plot) = cluster_exp{i_win}.cluster_size(i_cluster);
            end
        end
        
        %Plot topo
        cfg = [];
        cfg.layout    = 'CTF275.lay';
        cfg.colorbar  = 'yes';
        cfg.comment   = 'no';
        cfg.zlim = [-1 1];
        
        if n_cluster2plot ~= 0
            cfg.highlight        = 'on';
            cfg.highlightchannel = dat.label(cluster_sensorindices);
            cfg.highlightsymbol  = '.';
            cfg.highlightcolor   = [0.9 0.9 0.9];
            cfg.highlightsize    = 30;
        end
        
        h=figure;
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        fontsize = 15;
        
        ft_topoplotER(cfg, dat);
        
        %Determine cluster size
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
        
        %place title
        title({['Correlation Kprime & Rho-value (Prediction Performance) for behavioral p34 x p*34 ' ...
            'interaction effect across subjects'],...
            ['Kprime TW = [' w1 ' ms - ' w2 ' ms]; p-val < ' ...
            num2str(param_NASTD_MEG.plot.pval) ' (cluster-corrected); cluster size = ' cluster_size_text 'sensors']}, ...
            'FontSize', fontsize);
        
        %Save figure
        if param_NASTD_MEG.plot.save   == 1
            path_figs = [paths_NASTD_MEG.Current_outputfig 'Topo_ClusterCorr/'];
            mkdir(path_figs)
            
            filename     = ['Group_TopoClusterCorr_KprimeRho_' w1 '_' w2 'ms.png'];
            figfile      = [path_figs filename];
            saveas(gcf, [figfile], 'png'); %save png version
            delete(h);
        end    
    end
end

end