function [clusters_orig, clusters_shuffle, shuffleMaxStat] = NASTD_MEG_Pred_ClusterCorrect_Subs(dv_orig, iv_orig, pval_cluster, nReps_cluster)
% [clusters_orig, clusters_shuffle, shuffleMaxStat] = cluster_test_corr(dv, iv, p_thresh, fn_handle, nReps)
% 
% dv is a cell array of format
% dv{subject #}{sensor #}(trial #, 1)
%
% iv is a cell array of format
% iv{subject #}(trial #, iv#)
%
% fn_handle is a function handle to a function for analyzing the across-trial 
% relationship b/t dv and iv for each subject.
%
% The function referred to by fn_handle must be of the form
%
% [test_stat, p_val] = fun(dv, iv)
%
% where test_stat is the relevant summary statistic for the dv/iv relationship
% (e.g. regression beta or Pearson correlation z-score) and p_val is its 
% corresponding p-value.
%
% Shuffling takes place across trials within each individual subject data.
% Then the summary correlation statistic test_stat is subjected to an 
% across-subject one-sample t-test. The t value and p-value from this
% t-test are then used for cluster based analysis. Clusters are defined as
% sets of neighboring sensors such that every sensor in the cluster has a
% p-value less than p_thresh.


%% 1) Analysis of original data set 

nSens = length(dv_orig);
%1.1) Compute sensor-wise statistical effect
for i_sensor = 1:nSens

	stats_orig = regstats(dv_orig{i_sensor}, iv_orig, 'linear', 'tstat'); %linear regression of p33 MEG data on p*34       
    
    topo_p(i_sensor)    = stats_orig.tstat.pval(2);
    topo_stat(i_sensor) = stats_orig.tstat.t(2);
    
    stats_orig = [];

end

%1.2) Find sensor-clusters in the original data topo
clusters_orig = find_clusters(topo_stat, topo_p, pval_cluster);


%% 2) Analysis of repeatedly shuffled data set

for i_rep = 1:nReps_cluster    
    
    %2.1 Construct shuffle permutation for current iteration
    i_shuffledTrials = [];
	nTrials = length(dv_orig{1});% dv{sensor #}(trial #, 1)
	i_shuffledTrials = randperm(nTrials);

    %2.2 Apply shuffled trial index to  original data
    for i_sensor = 1:nSens
        dv_shuff{i_sensor} = dv_orig{i_sensor}(i_shuffledTrials);
        
        %2.3 Compute sensor-wise statistical effect
        stats_shuff = regstats(dv_shuff{i_sensor}, iv_orig, 'linear', 'tstat'); %linear regression of trial-shuffled p33 MEG data on p*34
        
        topo_p_shuffle(i_sensor)    = stats_shuff.tstat.pval(2);
        topo_stat_shuffle(i_sensor) = stats_shuff.tstat.t(2);
        
        stats_shuff = [];    
    end

    %2.4) Find clusters in the shuffled data topo and get the max cluster stat
    clusters_shuffle{i_rep} = find_clusters(topo_stat_shuffle, topo_p_shuffle, pval_cluster);
    shuffleMaxStat(i_rep)   = clusters_shuffle{i_rep}.maxStatSumAbs;

end


%% 3) Calculate p-values for each cluster in the original data set 
%  on the basis of the estimated null distribution from the shuffled data set
for i = 1:clusters_orig.nClusters
    pval = sum(shuffleMaxStat >= abs( clusters_orig.cluster_statSum(i) )) / nReps_cluster;
    clusters_orig.cluster_pval(i) = pval;
end