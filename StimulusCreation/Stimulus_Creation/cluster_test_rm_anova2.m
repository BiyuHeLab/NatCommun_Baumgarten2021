function [clusters_orig, clusters_shuffle, shuffleMaxStat] = cluster_test_rm_anova2(dv_orig, f_sub, f1, f2, f_names, f_analysis, p_thresh, nReps)
% [clusters_orig, clusters_shuffle, shuffleMaxStat] = cluster_test(dv, p_thresh, analysis, nReps)
% 
% format of dv
% 1 x nSensors cell array
% each cell holds a matrix of size (nSubs, nDVs) where nDVs correspond to the number of within-subject experimental conditions
%
% analysis is a string specifying the Matlab commands for conducting the
% statistical test on dv{i_sensor}. The depdendent variable in these commands 
% must be referred to by the variable name "dv{i_sensor}". The resulting
% p-value must be named "p_val" and the resulting test statistic must be
% named "test_stat".
%
% examples of constructing the analysis string
% ---------------------------------------------
% for repeated measures ANOVA:
%
%         analysis = [];
%         analysis = [analysis, '[p_val, table] = anova_rm(dv{i_sensor}, ''off'');  '];
%         analysis = [analysis, 'p_val = p_val(1);  '];
%         analysis = [analysis, 'test_stat = table{2, 5};  '];


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
nDVs     = size(dv_orig{1}, 2);

% compute topo of the statistical effect across subjects
for i_sensor = 1:nSensors

    stats = rm_anova2(dv_orig{i_sensor}(:), f_sub(:), f1(:), f2(:), f_names);
    p_val  = stats{f_analysis+1, 6};
    test_stat = stats{f_analysis+1, 5};

    topo_p(i_sensor)    = p_val;
    topo_stat(i_sensor) = test_stat;

end

% find clusters in the original data topo
clusters_orig = find_clusters(topo_stat, topo_p, p_thresh);



%% analysis of repeatedly shuffled data set

% shuffle

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
%             f1_shuffle(i_sub, :) = f1(i_sub, ind);
%             f2_shuffle(i_sub, :) = f2(i_sub, ind);
        end

%         stats = rm_anova2(dv_shuffle{i_sensor}(:), f_sub(:), f1_shuffle(:),f2_shuffle(:), f_names);
        stats = rm_anova2(dv_shuffle{i_sensor}(:), f_sub(:), f1(:),f2(:), f_names);
        p_val  = stats{f_analysis+1, 6};
        test_stat = stats{f_analysis+1, 5};

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