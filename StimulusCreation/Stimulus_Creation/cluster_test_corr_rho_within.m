function [clusters_orig, clusters_shuffle, shuffleMaxStat] = cluster_test_corr_rho_within(dv_orig, iv_orig, p_thresh, nReps)
% [clusters_orig, clusters_shuffle, shuffleMaxStat] = cluster_test_corr_rho_within(dv, iv, p_thresh, nReps)
% 
% Cluster based permutation correct for within-subject Spearman's rho
% correlation. Rho is computed across trials for each subject. The
% significance of rho being non-zero is tested across subjects by using an
% adjustment of the Fisher r-to-z transform for Spearman's rho: 
%
% z_rho = sqrt( (n-3)/1.06 ) * .5 * (log(1+rho) - log(1-rho))
% 
% 
% dv is a cell array of format
% dv{subject #}{sensor #}(trial #, 1)
%
% iv is a cell array of format
% iv{subject #}{sensor #}(trial #, 1)
%
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
        rho = corr(dv_orig{i_sub}{i_sensor}, iv_orig{i_sub}{i_sensor}, 'type', 'Spearman');
        
        n = length(dv_orig{i_sub}{i_sensor});
        test_stat(i_sub) = rho2z(rho, n);
    end
    
    [h, p_val, ci, stats] = ttest(test_stat);
    
    topo_p(i_sensor)    = p_val;
    topo_stat(i_sensor) = stats.tstat;

end

% find clusters in the original data topo
clusters_orig = find_clusters(topo_stat, topo_p, p_thresh);



%% analysis of repeatedly shuffled data set


% h = waitbar(0, 'Shuffling...');
for i_rep = 1:nReps

    waitbar(i_rep / nReps);
    
    % construct shuffle permutation for current iteration
    shuffle_ind = [];
    for i_sub = 1:nSubs
        % dv{subject #}{sensor #}(trial #, 1)
        nTrials = length(dv_orig{i_sub}{1});
        shuffle_ind{i_sub} = randperm(nTrials);
    end


    % 
    for i_sensor = 1:nSensors

        % apply the shuffling to the original data
        for i_sub = 1:nSubs
            ind = shuffle_ind{i_sub};
            dv_shuffle{i_sub}{i_sensor} = dv_orig{i_sub}{i_sensor}(ind);

%             test_stat(i_sub) = fn_handle(dv_shuffle{i_sub}{i_sensor}, iv_orig{i_sub}{i_sensor});
            
            rho = corr(dv_shuffle{i_sub}{i_sensor}, iv_orig{i_sub}{i_sensor}, 'type', 'Spearman');

            n = length(dv_shuffle{i_sub}{i_sensor});
            test_stat(i_sub) = rho2z(rho, n);            
        end

        [h, p_val, ci, stats] = ttest(test_stat);

        topo_p_shuffle(i_sensor)    = p_val;
        topo_stat_shuffle(i_sensor) = stats.tstat;
    
    end

    % find clusters in the shuffled data topo and get the max cluster stat
    clusters_shuffle{i_rep} = find_clusters(topo_stat_shuffle, topo_p_shuffle, p_thresh);
    shuffleMaxStat(i_rep)   = clusters_shuffle{i_rep}.maxStatSumAbs;

end

% delete(h);

%% calculate p-values for each cluster in the original data set 
%  on the basis of the estimated null distribution from the shuffled data set

for i = 1:clusters_orig.nClusters
    pval = sum(shuffleMaxStat > abs( clusters_orig.cluster_statSum(i) )) / nReps;
    clusters_orig.cluster_pval(i) = pval;
end