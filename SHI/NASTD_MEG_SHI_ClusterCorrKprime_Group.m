function NASTD_MEG_SHI_ClusterCorrKprime_Group...
    (subs, tonedur_text, ...
    nReps_GAvgNullDist, pval_clusterdef, ...
    plot_figs, pval_clusterplot, save_figs, ...
    paths_NASTD_MEG)

%Aims:
%1) Create group-level k-prime null-distribution from single-subject shuffled k-prime values
%2) Compute group-level statistical (cluster-corrected) k-prime effects (exp vs. shuff)
%3) Plot group-level statistical k-prime effects

%% 0) Specify vars, paths, and setup fieldtrip
addpath(genpath(paths_NASTD_MEG.ScriptsDir));

%Tone duration condition for data load-in
if strcmp(tonedur_text,'0.15')
    tonedur_title = '0.15sTD';
elseif  strcmp(tonedur_text,'0.3')
    tonedur_title = '0.3sTD';
elseif strcmp(tonedur_text,'0.6')
    tonedur_title = '0.6sTD';
end

path_outputdata = [paths_NASTD_MEG.Current_outputdata 'ExpK/' tonedur_title '/Group/'];
mkdir(path_outputdata)

%% 1) Define analysis parameters
win_size    = 30;
win_overlap = 0;

samplingFreq = 600;
toneDur_inSecs  = str2num(tonedur_text)
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

%% 2) Restructure single-subject Kprime-values in matrix containing all subs
for i_sub = 1:length(subs)
    path_inputdata = [paths_NASTD_MEG.Current_outputdata 'ExpK/' tonedur_title '/' subs{i_sub} '/'];
    load([path_inputdata subs{i_sub} '_ExpvsShuffKprime_' tonedur_title '.mat' ]);
    
    for i_win = 1:nWindows
        %Create summary file with experimental k-prime vals for all subs
        KprimeEXP_AllSub{i_win}(i_sub, :) = Kprime.Exp.avgFolds.Kprime.Avg{i_win}; %Exp k-prime vals (averaged over folds)
        %Create summary file with shuffled K-prime vals for all subs
        KprimeSHUFFLED_AllSub{i_win}{i_sub} = Kprime.Shuff.avgFolds.NullDist{i_win}'; 
        %shuffled k-prime vals (for each rep, averaged across folds)
    end
    clear Kprime
end

%% 3) Create shuffled data null distribution
nSensors       = size(KprimeEXP_AllSub{1},2);
nWithinSubReps = size(KprimeSHUFFLED_AllSub{1}{1},1);

%Average shuffled k-prime values across subjects
for i_rep = 1:nReps_GAvgNullDist
    %Repeated random sampling (with replacement) from the permuted
    %k-prime values (100 reps per fold, each rep averaged across folds
    %= 100 reps total per subject) of each subject/time window/sensor.
    %Calculation of across-subject average for 'nReps_GAvgNullDist'time
    %yielding a distribution of mean k-prime values under the
    %hypothesis that k-prime = 0.
    
    for i_sub = 1:length(subs)
        rand_ind = ceil(rand*nWithinSubReps); %Create random index between 1-100 for each rep and subject
        
        for i_win = 1:nWindows
            %Read out one of the 100 all-sensor arrays of shuffled kprime values randomly determined
            %by random index for each subject, per time window and repeptition
            KprimeSHUFFLED_RandSelection_AllSub{i_rep}{i_win}(i_sub, :) = KprimeSHUFFLED_AllSub{i_win}{i_sub}(rand_ind, :);
            %Output: Take one random shuffled k-prime per sensor array per subject for each time window and repetition
        end
    end
    
    for i_win = 1:nWindows
        %Average the randomly selected k-prime values across subjects for each
        %rep, to create Gavg null distribution with 1000 Gavg entries
        KprimeSHUFFLED_RandSelection_GAvg{i_win}(i_rep, :) = mean(KprimeSHUFFLED_RandSelection_AllSub{i_rep}{i_win});
        %Output: Across-subject average of shuffled k-prime values per
        %sensor for each repetition and time window
    end
end

%Calculate p-values according to shuffled k-prime null distribution
for i_win = 1:nWindows
    %Average experimental k-values across subjects for each sensor and time window
    KprimeEXP_GAvg{i_win} = mean(KprimeEXP_AllSub{i_win});
    
    for i_sensor = 1:nSensors
        %Determine the number of sensors where shuffled k-prime values
        %from the 1000 across-subject averages are larger than or equal
        %to the experimental k-prime
        Kprime_EXPvsSHUFFLED_GAvg_pval{i_win}(i_sensor) = ...
            sum(KprimeSHUFFLED_RandSelection_GAvg{i_win}(:, i_sensor) ...
            >= KprimeEXP_GAvg{i_win}(i_sensor)) / nReps_GAvgNullDist;
        
        %For each repetition, determine the number of sensors where
        %shuffled k-prime values of that repetition are larger than
        %shuffled k-prime values from the 1000 across-subject averages
        for i_rep = 1:nReps_GAvgNullDist
            KprimeSHUFFLED_PermvsAvgPerm_pval{i_win}(i_rep, i_sensor) = ...
                sum(KprimeSHUFFLED_RandSelection_GAvg{i_win}(:, i_sensor)...
                >= KprimeSHUFFLED_RandSelection_GAvg{i_win}(i_rep, i_sensor))...
                / nReps_GAvgNullDist;
        end
    end
end

%% 4) save Output
savefile = [path_outputdata 'Group_ShuffKprimeNullDist_' tonedur_title '.mat'];
save(savefile, ...
    'KprimeEXP_GAvg', ...
    'KprimeSHUFFLED_RandSelection_AllSub', 'KprimeSHUFFLED_RandSelection_GAvg', ...
    'Kprime_EXPvsSHUFFLED_GAvg_pval', 'KprimeSHUFFLED_PermvsAvgPerm_pval');

%% 5) Compute Group-level cluster statistics
%Perform Cluster correction
nReps_GavgNullDist = size(KprimeSHUFFLED_RandSelection_GAvg{1},1);
%Across-subject average of shuffled k-prime values per
%sensor for each repetition and time window
nSubs = size(KprimeSHUFFLED_RandSelection_AllSub{1}{1},1);
%Take one random shuffled k-prime per sensor array per
%subject for each time window and repetition

data_cluster_filename = ['Group_KprimeClusterCorr_' tonedur_title ...
    '_Clusterp' num2str(pval_clusterdef) '_nReps' num2str(nReps_GavgNullDist) ...
    '_nSubs' num2str(nSubs) '.mat'];

for i_win = 1:nWindows
    
    %find clusters in the original data topo
    topo_GAvgEXPKval = KprimeEXP_GAvg{i_win};
    %Experimental k-prime-values averaged across subjects for each sensor and time window
    topo_GAvgEXPKval_pval = Kprime_EXPvsSHUFFLED_GAvg_pval{i_win};
    %Experimental vs shuffled null distribution p-values for each sensor and time window
    
    clusters_EXP{i_win} = find_clusters(topo_GAvgEXPKval, topo_GAvgEXPKval_pval, pval_clusterdef);
    
    % find clusters in the shuffled data topo and get the max cluster stat to compare against the Exp cluster
    for i_rep = 1:nReps_GavgNullDist
        topo_SHUFFLEDKval_Gavg        = KprimeSHUFFLED_RandSelection_GAvg{i_win}(i_rep, :);
        topo_SHUFFLEDKval_Gavg_pval   = KprimeSHUFFLED_PermvsAvgPerm_pval{i_win}(i_rep, :);
        
        clusters_SHUFFLED  = find_clusters(topo_SHUFFLEDKval_Gavg, topo_SHUFFLEDKval_Gavg_pval, pval_clusterdef);
        MaxStat_SHUFFLED{i_win}(i_rep) = clusters_SHUFFLED.maxStatSumAbs;
    end
    
    % calculate p-values for clusters in original data using clusters computed with shuffled data
    for i = 1:clusters_EXP{i_win}.nClusters %for each cluster found inthe original data
        pval = sum(MaxStat_SHUFFLED{i_win} >= abs( clusters_EXP{i_win}.cluster_statSum(i) )) / nReps_GavgNullDist;
        %For each cluster, compute number of max cluster statistic in shuffled data that are
        %bigger than the max cluster statistic of the experimental
        %condition divided by number of repritions
        clusters_EXP{i_win}.cluster_pval(i) = pval;
    end
    
end

%% 6) Save Group-level cluster output
save([path_outputdata, data_cluster_filename], 'clusters_EXP', 'MaxStat_SHUFFLED', 'nReps_GavgNullDist', '-v7.3');

%% 7) Plot topo with Group-level cluster-corrected k-prime values
if plot_figs == 1
    load([paths_NASTD_MEG.ScriptsDir 'MEG_sensor_setup_272/label272.mat']); %file with CTF sensor labels for 272 sensors, called 'label'
    
    %Determine z scaling
    absmax_plot = -Inf;
    absmin_plot = +Inf;
    
    for i_win = 1:nWindows
        GAvgCluster_topoplot = clusters_EXP{i_win}.inputs.topo_stat;
        
        meanKprime(i_win) = mean(clusters_EXP{i_win}.inputs.topo_stat);
        stdKprime(i_win) = std(clusters_EXP{i_win}.inputs.topo_stat);
        
        if max(abs(GAvgCluster_topoplot)) > absmax_plot
            absmax_plot = floor(max(abs(GAvgCluster_topoplot)));
        end
        if min(abs(GAvgCluster_topoplot)) < absmin_plot
            absmin_plot = floor(min(abs(GAvgCluster_topoplot)));
        end
    end
    
    %Create topo struct with group-avg Exp k-prime vals per sensor for current time window
    
    for i_win = 1:nWindows
        w1 = num2str( 1000 * windows(i_win, 1) / samplingFreq , 3 );
        w2 = num2str( 1000 * windows(i_win, 2) / samplingFreq , 3 );
        
        nSensors = length(clusters_EXP{1}.topo_cluster);
        GAvgCluster_topoplot = clusters_EXP{i_win}.inputs.topo_stat; %Group-avg Exp k-prime values per sensor
        
        clear dat
        dat.dimord = 'chan_time';
        dat.label  = label;
        dat.time   = 0;
        dat.avg = GAvgCluster_topoplot';
        
        %Determine if cluster is significant and thus plotted
        n_cluster2plot = 0;
        cluster_sensorindices = [];
        
        for i_cluster = 1:clusters_EXP{i_win}.nClusters%for each found cluster
            cluster_pval = clusters_EXP{i_win}.cluster_pval(i_cluster); %Read out cluster-specific p-val
            
            if cluster_pval < pval_clusterplot %if p-val is below plotting threshold, add cluster to plot list
                n_cluster2plot = n_cluster2plot + 1;
                cluster_sensorindices = [cluster_sensorindices clusters_EXP{i_win}.cluster_sensors{i_cluster}];
                %Read out sensor indices for significant clusters (appended for multiple clusters)
                cluster_size(n_cluster2plot) = clusters_EXP{i_win}.cluster_size(i_cluster);
                %Read out cluster size (number of sensors in cluster)
            end
        end
        
        %3.4 Plot topo
        cfg = [];
        cfg.layout    = 'CTF275.lay';
        cfg.colorbar  = 'yes';
        cfg.comment   = 'no';
        cfg.zlim = [absmin_plot, absmax_plot];
        cfg.zlim = [4, 8];
        
        cfg.highlight        = 'on';
        cfg.highlightchannel = dat.label(cluster_sensorindices);
        cfg.highlightsymbol  = '.';
        cfg.highlightcolor   = [0.9 0.9 0.9];
        cfg.highlightsize    = 30;
        
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
        
        h = figure;
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        fontsize = 15;
        
        ft_topoplotER(cfg, dat);%
        title({['Group-level Sign. Kprime (n = ' num2str(nSubs) '; Cluster-corrected)'],...
            [tonedur_title '; ERF TW = [' w1 ' ms - ' w2 ' ms]; '], ...
            [num2str(nReps_GAvgNullDist) 'Reps; Cluster-pval: ' num2str(pval_clusterplot) '; cluster sizes = ' cluster_size_text]}, ...
            'FontSize', fontsize);
        
        if save_figs    == 1
            path_figs = [paths_NASTD_MEG.Current_outputfig 'ExpK/' tonedur_title '/Group/'];
            mkdir(path_figs)
            
            filename     = ['Group_SignKprime_ClusterCorr_' tonedur_title '_' w1 '-' w2 '.png'];
            figfile      = [path_figs filename];
            saveas(gcf, [figfile], 'png'); %save png version
            delete(h);
        end
        
    end
end

end
