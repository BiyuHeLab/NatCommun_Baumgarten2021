function clusterSOI = NASTD_MEG_ERF_GetSensorIndexfromCluster(tonedur_text, param, paths_NASTD_MEG)

%Script reads out sensor indices for particular predictive processing sensor cluster

%% 1) Determine Input data (predictive processing luster data per tone duration)
path_inputdata = ...
    [paths_NASTD_MEG.Analysis.Prediction ...
        'Data/GAvgStatistics_ClusterCorr/not_baseline_corrected/'...
        num2str(str2num(tonedur_text)*1000) 'ms_TD/']; %old path where GAvg statistics are saved
%     [paths_NASTD_MEG.Analysis.Prediction 'PredEffect/' tonedur_text 'sTD/Group/']; %new path where GAvg statistics are saved

load([path_inputdata 'ClusterPred_' num2str(str2num(tonedur_text)*1000) ...
    'msTD_p' num2str(param.SOI.pval_load) '_nReps' num2str(param.SOI.nPerm) ...
    '_nSubs' num2str(length(param.subs)) '.mat']);

% load([path_inputdata 'Group_PredEffect_' tonedur_text ...
%     'sTD_ClusterCorrect_p' num2str(param.cluster.pval_load) '_nReps' num2str(param.cluster.nPerm) '.mat']);

%% 2) Read out Sensors
clusterSOI = [];
numCluster = 0;
SensperCluster =[];

for i_win = param.SOI.windows %Loop across timewin (lowest common window number)    
    if clusters_Exp{i_win}.nClusters > 0 %in case there are sign clusters present
        i_signClusters = ...
            find(clusters_Exp{i_win}.cluster_pval < param.SOI.pval_select); %find which clusters are considered sign.
        
        if ~isempty(i_signClusters) %if there are any sign. clusters
            for i_clusters = 1:length(i_signClusters)
                numCluster = numCluster+1;
                clusterSOI{i_win}.i_SensperCluster{i_clusters} = ...
                    clusters_Exp{i_win}.cluster_sensors{i_signClusters(i_clusters)};
                %Read out sensors of sign. cluster
                SensperCluster = ...
                    [SensperCluster, length(clusterSOI{i_win}.i_SensperCluster{i_clusters})];
                
                clusterSOI{i_win}.clusterSign{i_clusters} = ...
                    clusters_Exp{i_win}.cluster_statSum(i_signClusters(i_clusters)) > 1;
                %Check if cluster effect is positive or negative
            end
        else
            clusterSOI{i_win}.i_SensperCluster = [];
            clusterSOI{i_win}.clusterSign = [];
        end
    end
    
    %Place all sensors of all clusters in common vector
    clusterSOI{i_win}.i_allSens = [];
    if ~isempty(clusterSOI{i_win}.i_SensperCluster) %if there are any sign. clusters
        for i_cluster = 1:length(clusterSOI{i_win}.i_SensperCluster) %loop across clusters
            clusterSOI{i_win}.i_allSens = ...
                [clusterSOI{i_win}.i_allSens clusterSOI{i_win}.i_SensperCluster{i_cluster}];
        end
    end
    
end
disp(['Found: ' num2str(numCluster) ' sign. Clusters in total with ' num2str(round(mean(SensperCluster))) ' Sens on average.'])

end