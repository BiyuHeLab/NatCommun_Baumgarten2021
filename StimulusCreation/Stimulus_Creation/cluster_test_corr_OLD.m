function [clusters_orig, clusters_shuffle, shuffleMaxStat] = cluster_test_corr(dv_orig, iv_orig, p_thresh, fn_handle, nReps)
% [clusters_orig, clusters_shuffle, shuffleMaxStat] = cluster_test_corr(dv, iv, p_thresh, fn_handle, nReps)
% 
% format of dv
% 1 x nSensors cell array
% each cell holds a matrix of size (nSubs, nDVs) where nDVs correspond to the number of within-subject experimental conditions
%
% fn_handle is a function handle to a function for analyzing DV and IV.
% The function must be of the form
%
% [test_stat, p_val] = fn_handle(dv, iv)
%
% where test_stat is the relevant test statistic for cluster-based analysis
% and p_val is its corresponding p-value.



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
nSubs    = size(dv_orig{1}, 1);

% compute topo of the statistical effect across subjects
for i_sensor = 1:nSensors

    
%     [pval, table] = anova_rm(dv{i_sensor}, 'off');
%     pval = pval(1);
%     stat = table{2, 5};

    dv = dv_orig;
    iv = iv_orig;
    
    [test_stat, p_val] = fn_handle(dv, iv);
    
    topo_p(i_sensor)    = p_val;
    topo_stat(i_sensor) = test_stat;

end

% find clusters in the original data topo
clusters_orig = find_clusters(topo_stat, topo_p, p_thresh);



%% analysis of repeatedly shuffled data set

% shuffle
nSubs = size(dv{1}, 1);
nDVs  = size(dv{1}, 2);

h = waitbar(0, 'Shuffling...');
for i_rep = 1:nReps

    waitbar(i_rep / nReps);
    
    % construct shuffle permutation for current iteration
    shuffle_ind = [];
    for i_sub = 1:nSubs
        shuffle_ind(i_sub, :) = randperm(nDVs);
    end


    % 
    for i_sensor = 1:nSensors

        % apply the shuffling to the original data
        for i_sub = 1:nSubs
            ind = shuffle_ind(i_sub, :);
            dv_shuffle{i_sensor}(i_sub, :) = dv_orig{i_sensor}(i_sub, ind);
        end

        % compute stat topo for shuffled set
    %     [pval, table] = anova_rm(dv{i_sensor}, 'off');
    %     pval = pval(1);
    %     stat = table{2, 5};

        dv = dv_shuffle;
        [test_stat, p_val] = fn_handle(dv, iv);

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
    pval = sum(shuffleMaxStat > abs( clusters_orig.cluster_statSum(i) )) / nReps;
    clusters_orig.cluster_pval(i) = pval;
end