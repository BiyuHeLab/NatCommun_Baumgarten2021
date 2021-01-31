function [clusters_orig, clusters_shuffle, shuffleMaxStat] = NASTD_MEG_Pred_ClusterCorrect_Group(dv_orig, iv_orig, p_thresh, nReps)
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

%% analysis of original data set

nSubs    = length(dv_orig);
nSensors = length(dv_orig{1});

% compute topo of the statistical effect across subjects
for i_sensor = 1:nSensors

    for i_sub = 1:nSubs
        stats1 = regstats(dv_orig{i_sub}{i_sensor}, iv_orig{i_sub}, 'linear', 'tstat');
        test_stat(i_sub) = stats1.tstat.beta(2); %Beta value of linear regression as summary stat
        stats1 = [];
    end
    
    [h, p_val, ci, stats] = ttest(test_stat);
    
    topo_p(i_sensor)    = p_val;
    topo_stat(i_sensor) = stats.tstat;

end

% find clusters in the original data topo
clusters_orig = find_clusters(topo_stat, topo_p, p_thresh);

%% analysis of repeatedly shuffled data set
for i_rep = 1:nReps
    disp(['rep#: ' num2str(i_rep) '/' num2str(nReps)])
    % construct shuffle permutation for current iteration
    shuffle_ind = [];
    for i_sub = 1:nSubs
        nTrials = length(dv_orig{i_sub}{1});
        shuffle_ind{i_sub} = randperm(nTrials);
    end
    
    for i_sensor = 1:nSensors

        % apply the shuffling to the original data
        for i_sub = 1:nSubs
            ind = shuffle_ind{i_sub};
            dv_shuffle{i_sub}{i_sensor} = dv_orig{i_sub}{i_sensor}(ind);
            stats1 = regstats(dv_shuffle{i_sub}{i_sensor}, iv_orig{i_sub}, 'linear', 'tstat');
            test_stat(i_sub) = stats1.tstat.t(2);
        end

        [h, p_val, ci, stats] = ttest(test_stat);

        topo_p_shuffle(i_sensor)    = p_val;
        topo_stat_shuffle(i_sensor) = stats.tstat;
    
    end

    % find clusters in the shuffled data topo and get the max cluster stat
    clusters_shuffle{i_rep} = find_clusters(topo_stat_shuffle, topo_p_shuffle, p_thresh);
    shuffleMaxStat(i_rep)   = clusters_shuffle{i_rep}.maxStatSumAbs;

end

%% calculate p-values for each cluster in the original data set 
%  on the basis of the estimated null distribution from the shuffled data set
for i = 1:clusters_orig.nClusters
    pval = sum(shuffleMaxStat >= abs( clusters_orig.cluster_statSum(i) )) / nReps;
    clusters_orig.cluster_pval(i) = pval;
end