function NASTD_MEG_SHI_ClusterCorrKprime_Group_perBeta...
    (subs, betalevel_input, tonedur_text, ...
    nReps_GAvgNullDist, pval_clusterdef, ...
    pval_clusterplot, save_figs, ...
    paths_NASTD_MEG)

%Aims:
%1) Create group-level k-prime null-distribution from single-subject shuffled k-prime values
%2) Compute group-level statistical (cluster-corrected) k-prime effects (exp vs. shuff)
%3) Plot group-level statistical k-prime effects

%% 0) specify paths and set up vars
addpath(genpath(paths_NASTD_MEG.ScriptsDir));

nSubs = length(subs);

%Data input/output
path_inputdata = [paths_NASTD_MEG.Current_outputdata 'ExpK/SeqBeta_' betalevel_input '/'];
path_outputdata = [paths_NASTD_MEG.Current_outputdata 'ExpK/SeqBeta_' betalevel_input '/Group/'];
path_fig =  [paths_NASTD_MEG.Current_outputfig 'ExpK/SeqBeta_' betalevel_input '/Group/'];
mkdir(path_outputdata);

%% 1) Define analysis parameters
win_size    = 30;
win_overlap = 0;

samplingFreq = 600;
toneDur_inSecs  = str2num(tonedur_text{1});
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

%% 2) Restructure single-subject k-values in matrix containing all subs
for i_sub = 1:length(subs)    
              
    load([path_inputdata subs{i_sub} '_ExpvsShuffKprime_SeqBeta' betalevel_input '.mat' ]);

    for i_win = 1:nWindows        
        %Create summary file with experimental K-prime vals for all subs
        %(subs = rows, sensors = columns)
        KprimeEXP_AllSub{i_win}(i_sub, :) = ...
            Kprime.Exp.avgFolds.Kprime.Avg{i_win}; %Exp K-prime vals (averaged over folds)
        
        %Create summary file with shuffled K-prime vals for all subs
        %subs = cells, rows = repetitions, columns = sensors)
        KprimeSHUFFLED_AllSub{i_win}{i_sub} = ...
            Kprime.Shuff.avgFolds.NullDist{i_win}'; %shuffled K-prime vals (for each rep, averaged across folds)
    end
        clear Kprime
end

%% 3) Create shuffled data null distribution
nSensors       = size(KprimeEXP_AllSub{1},2);
nWithinSubReps = size(KprimeSHUFFLED_AllSub{1}{1},1); 
    
%3.1 Average shuffled k-prime values across subjects
for i_rep = 1:nReps_GAvgNullDist
%Note: Repeated random sampling (with replacement) from the permuted
%k-prime values (100 reps per fold, each rep averaged across folds 
%= 100 reps total per subject) of each subject/time window/sensor. 
%Calculation of across-subject average for x time (Maniscalco et al. 2018 = 
%1000 times) yielding a distribution of mean k-prime values under the 
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

%3.2 Calculate p-values according to shuffled null distribution
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

%3.3 save variables
savefile = [path_outputdata 'Group_ShuffKprimeNullDist_SeqBeta' betalevel_input '.mat'];

save(savefile, ...
    'KprimeEXP_GAvg', ...
    'KprimeSHUFFLED_RandSelection_AllSub', 'KprimeSHUFFLED_RandSelection_GAvg', ...
    'Kprime_EXPvsSHUFFLED_GAvg_pval', 'KprimeSHUFFLED_PermvsAvgPerm_pval');

%% 4) Perform Cluster correction
nReps_GavgNullDist = size(KprimeSHUFFLED_RandSelection_GAvg{1},1);
            %Across-subject average of shuffled k-prime values per
            %sensor for each repetition and time window
nSubs = size(KprimeSHUFFLED_RandSelection_AllSub{1}{1},1);
            %Take one random shuffled k-prime per sensor array per 
            %subject for each time window and repetition

for i_win = 1:nWindows
    
    % find clusters in the original data topo
    topo_GAvgEXPKval = KprimeEXP_GAvg{i_win}; %Experimental k-prime-values averaged across subjects for each sensor and time window
    topo_GAvgEXPKval_pval = Kprime_EXPvsSHUFFLED_GAvg_pval{i_win}; %Experimental vs shuffled null distribution p-values for each sensor and time window
    
    clusters_EXP{i_win} = ...
        find_clusters(topo_GAvgEXPKval, topo_GAvgEXPKval_pval, pval_clusterdef); 
        
    % find clusters in the shuffled data topo and get the max cluster stat
    % to compare against the Exp cluster
    for i_rep = 1:nReps_GavgNullDist
        topo_SHUFFLEDKval_Gavg      = KprimeSHUFFLED_RandSelection_GAvg{i_win}(i_rep, :);
        topo_SHUFFLEDKval_Gavg_pval = KprimeSHUFFLED_PermvsAvgPerm_pval{i_win}(i_rep, :);
        
        clusters_SHUFFLED  = find_clusters(topo_SHUFFLEDKval_Gavg, topo_SHUFFLEDKval_Gavg_pval,pval_clusterdef );
        MaxStat_SHUFFLED{i_win}(i_rep) = clusters_SHUFFLED.maxStatSumAbs;
    end

    % calculate p-values for clusters in original data using clusters computed with shuffled data
    for i = 1:clusters_EXP{i_win}.nClusters %for each cluster found inthe original data
        pval = sum(MaxStat_SHUFFLED{i_win} >= abs( clusters_EXP{i_win}.cluster_statSum(i) )) / nReps_GavgNullDist;
            %For each cluster, compute number of max cluster statistic in shuffled data that are
            %bigger than the max cluster statistic of the experimental
            %conditionm divided by number of repritions
        clusters_EXP{i_win}.cluster_pval(i) = pval;
    end    
    
end

%% 4) Save cluster output
data_cluster_filename = ['Group_KprimeClusterCorr_SeqBeta' betalevel_input ...
    '_Clusterp' num2str(pval_clusterdef) '_nReps' num2str(nReps_GavgNullDist) ...
    '_nSubs' num2str(nSubs) '.mat'];

save([path_outputdata, data_cluster_filename], 'clusters_EXP', 'MaxStat_SHUFFLED', 'nReps_GavgNullDist', '-v7.3');

%% 5) Plot topo with cluster corrected across subject average k-prime values
load([paths_NASTD_MEG.ScriptsDir 'MEG_sensor_setup_272/label272.mat'] ); 
%file with CTF sensor labels for 272 sensors

%Determine z scaling
absmax_plot = -Inf;
absmin_plot = +Inf;

for i_win = 1:nWindows
    GAvgCluster_topoplot = clusters_EXP{i_win}.inputs.topo_stat;
    
    meanKprime(i_win) = mean(clusters_EXP{i_win}.inputs.topo_stat);
    varKprime(i_win) = var(clusters_EXP{i_win}.inputs.topo_stat);
    stdKprime(i_win) = std(clusters_EXP{i_win}.inputs.topo_stat);
    
    if max(abs(GAvgCluster_topoplot)) > absmax_plot
        absmax_plot = floor(max(abs(GAvgCluster_topoplot)));
    end
    if min(abs(GAvgCluster_topoplot)) < absmin_plot
        absmin_plot = floor(min(abs(GAvgCluster_topoplot)));
    end    
end
SDKprime_acrossTW = sqrt(mean(varKprime));
meanKrpime_acrossTW = mean(meanKprime);

%Set up titles and Topo struct
model_comp_title = 'Kprime (defined by minSumSquaredResiduals/ModelOrder for linear regression on tone history) per sensor';
title_text = ['Exp vs. Shuff (ClusterCorrected, p = ' num2str(pval_clusterplot) ')'];

clear dat
dat.dimord = 'chan_time';
dat.label  = label;
dat.time   = 0;

cfg = [];
cfg.layout    = 'CTF275.lay';
cfg.comment   = 'no';
cfg.colorbar  = 'yes';
% cfg.zlim      = [absmin_plot, absmax_plot]; 
cfg.zlim      = [4, 8];

%Plot all topos/TW in one figure
pos_fig = 0;
figure; set(gcf,'units','normalized','outerposition',[0 0 1 1])

for i_win = 1:nWindows      
    
    %Determine if cluster is significant and thus plotted
    n_cluster2plot = 0;
    cluster_sensorindices = [];

    for i_cluster = 1:clusters_EXP{i_win}.nClusters%for each found cluster     
        cluster_pval = clusters_EXP{i_win}.cluster_pval(i_cluster); %Read out cluster-specific p-val
        
        if cluster_pval < pval_clusterplot %if p-val is below plotting threshold, add cluster to plot list
            n_cluster2plot = n_cluster2plot + 1;%add 1 to plot list counter            
            cluster_sensorindices = [cluster_sensorindices clusters_EXP{i_win}.cluster_sensors{i_cluster}];
            %Read out sensor indices for significant clusters (appended for
            %multiple clusters)
            cluster_size(n_cluster2plot) = clusters_EXP{i_win}.cluster_size(i_cluster);
            %Read out cluster size (number of sensors in cluster)         
        end        
    end    
    
    %determine start & end points for respective window
    w1 = num2str(1000 * windows(i_win, 1) / samplingFreq , 3 );
    w2 = num2str(1000 * windows(i_win, 2) / samplingFreq , 3 );
    win_title = ['TW = [' w1 ' ms - ' w2 ' ms]'];
    
    cfg.colorbar  = 'no';
    
    colormap('parula')
    clim = [absmin_plot, absmax_plot];    
    clim = [4, 8];
   
    cfg.style           = 'both'; % 'straight_imsat'; %reduces file size
    cfg.marker          = 'on';
    cfg.markersymbol    = '.';
    cfg.markersize      = 6;
    
    cfg.highlight = 'on';
    cfg.highlightchannel = label(cluster_sensorindices);
    %Highlight those channels where kprime is significantly different
    %from shuffled null distribution
    %i.e., dertermined by the amount of shuffled k-prime values that are
    %equal or bigger to the experimental/original k-prime value and
    %divides this by the number of reps
    cfg.highlightcolor      = [0.9 0.9 0.9];
    cfg.highlightsymbol     = '.';
    cfg.highlightsize       = 20;
    
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
         
    dat.avg = clusters_EXP{i_win}.inputs.topo_stat'; %Group-avg Exp k-prime vals per sensor    
    
    pos_fig = pos_fig + 1;
    hold on;
    subplot(nWindows/3, 3, pos_fig);
    ft_topoplotER(cfg, dat);
    t = title({[win_title],['Cluster Size: ' cluster_size_text ' channels']});
    t.Position(2) = t.Position(2)-(t.Position(2)*0.2);
    
    if i_win == nWindows
        h = colorbar;
        h.Limits = clim;
        h.Position(1) = 0.91; %sets colorbar higher
        h.Position(2) = 0.225; %sets colorbar higher
        h.Position(4) = h.Position(4)*0.7; %makes colorbar shorter
        h.Ticks = [4 5 6 7 8];
        h.Label.String = 'Exp Kprime';       
        set(h,'FontSize',16);
    end
       
end
suptitle({['GAvg (n = ' num2str(nSubs) ') - SeqBeta: ' betalevel_input ], model_comp_title, title_text});

if save_figs
    mkdir(path_fig); %make new directory for output files
    filename     = ['Group_SignKprime_ClusterCorr_SeqBeta' betalevel_input '_AllTW.png'];
    figfile      = [path_fig filename];
    
    saveas(gcf, [figfile], 'png'); %save png version
    close
end

end