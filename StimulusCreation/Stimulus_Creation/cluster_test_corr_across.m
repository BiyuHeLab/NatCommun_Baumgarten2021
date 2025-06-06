function [clusters_orig, clusters_shuffle, shuffleMaxStat] = cluster_test_corr_across(dv_orig, p_thresh, corr_type, nReps)
% [clusters_orig, clusters_shuffle, shuffleMaxStat] = cluster_test_corr_across(dv, p_thresh, analysis, nReps)
% 
% Cluster based permutation correct for across-subject correlation.
%
% format of dv
% 1 x nSensors cell array
% each cell holds a matrix of size (nSubs, 2) where the two columns correspond 
% to the 2 variables being correlated across subjects.
%
% p_thresh is the p-value used as the threshold for defining sensor clusters.
%
% corr_type is the type of correlation being conducted (e.g. 'Pearson' or
% 'Spearman' or 'Kendall').
%
% nReps is the number of repetitions used for permutation shuffling.


% original data
% -------------
% for each sensor,
% - compute topo of the statistical effect across subjects
% across all sensors,
% - find clusters in the original data topo
%
% shuffled data
% -------------
% for each sensor,
% - shuffle experimental condition means randomly and independently for each subject
% across all sensors,
% - compute topo of the statistical effect across subjects for the shuffled data set
% - find and save the max(ClusterStat) for thus shuffled set
% - over many iterations, this gives the null distribution for ClusterStat


%% analysis of original data set

nSensors = length(dv_orig);

% compute topo of the statistical effect across subjects
for i_sensor = 1:nSensors

    
%     [pval, table] = anova_rm(dv{i_sensor}, 'off');
%     pval = pval(1);
%     stat = table{2, 5};

    dv = dv_orig;
    [test_stat, p_val] = corr(dv{i_sensor}(:,1), dv{i_sensor}(:,2), 'type', corr_type);

    topo_p(i_sensor)    = p_val;
    topo_stat(i_sensor) = test_stat;

end

% find clusters in the original data topo
clusters_orig = find_clusters(topo_stat, topo_p, p_thresh);



%% analysis of repeatedly shuffled data set

% shuffle
nSubs = size(dv{1}, 1);

h = waitbar(0, 'Shuffling...');
for i_rep = 1:nReps

    waitbar(i_rep / nReps);
    
    % construct shuffle permutation for current iteration
    shuffle_ind = randperm(nSubs);


    % 
    for i_sensor = 1:nSensors

        % apply the shuffling to the original data
        for i_sub = 1:nSubs
            dv_shuffle{i_sensor}(:, 1) = dv_orig{i_sensor}(:, 1);
            dv_shuffle{i_sensor}(:, 2) = dv_orig{i_sensor}(shuffle_ind, 2);
        end

        % compute stat topo for shuffled set
    %     [pval, table] = anova_rm(dv{i_sensor}, 'off');
    %     pval = pval(1);
    %     stat = table{2, 5};

        dv = dv_shuffle;
        [test_stat, p_val] = corr(dv{i_sensor}(:,1), dv{i_sensor}(:,2), 'type', corr_type);

        topo_p_shuffle(i_sensor)    = p_val;
        topo_stat_shuffle(i_sensor) = test_stat;
        
    end

    % find clusters in the shuffled data topo and get the max cluster stat
    clusters_shuffle{i_rep} = find_clusters(topo_stat_shuffle, topo_p_shuffle, p_thresh);
    shuffleMaxStat(i_rep) = clusters_shuffle{i_rep}.maxStatSumAbs;

end

delete(h);

%% calculate p-values for each cluster in the original data set 
%  on the basis of the estimated null distribution from the shuffled data set

for i = 1:clusters_orig.nClusters
    pval = sum(shuffleMaxStat > abs( clusters_orig.cluster_statSum(i)) )/ nReps;
    clusters_orig.cluster_pval(i) = pval;
end