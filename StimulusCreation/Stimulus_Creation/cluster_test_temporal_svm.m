function [clusters_orig, clusters_shuffle, shuffleMaxStatPos] = cluster_test_temporal_svm(labels, data, p_thresh, nReps, libsvm_settings)
% [clusters_orig, clusters_shuffle, shuffleMaxStatPos] = cluster_test_temporal_svm(group, data, p_thresh, nReps, libsvm_settings)
% 
% labels - 1 x nsubs cell array holding ntrials x 1 array of trial labels for each subject
% labels{nsubs}(ntrials, 1)
%
% data  - 1 x nsubs cell array holding 1 x nsamples cell array holding ntrials x nsensors data matrix 
% data{nsubs}{nsamples}(ntrials, nsensors)


% original data
% -------------
% for each sample,
% - compute timecourse of the statistical effect across subjects
% across all samples,
% - find clusters in the original data timecourse
%
% shuffled data
% -------------
% for each sample,
% - shuffle experimental condition means randomly and independently for each subject
% across all samples,
% - compute timecourse of the statistical effect across subjects for the shuffled data set
% - find and save the max(ClusterStat) for thus shuffled set
% - over many iterations, this gives the null distribution for ClusterStat

%%

nSubs    = length(labels);
nSamples = length(data{1});
nSensors = size(data{1}{1}, 2);

nfolds     = libsvm_settings.nfolds;
nreps_svm  = libsvm_settings.nreps_svm;
scale_data = libsvm_settings.scale_data;
libsvm_options_train = libsvm_settings.libsvm_options_train;
libsvm_options_test  = libsvm_settings.libsvm_options_test;



%% analysis of original data set

tic

% get SVM accuracy for each subject and sample
for i_sub = 1:nSubs

    labels_i = labels{i_sub};
    
    for i_sample = 1:nSamples
        
        data_i   = data{i_sub}{i_sample};
        
        svm_acc = libsvm_nfold(labels_i, data_i, nfolds, nreps_svm, libsvm_options_train, libsvm_options_test, scale_data);       
        
        svm_acc_time(i_sub, i_sample) = mean(svm_acc(:));
        
    end
    
end

% get statistical summary of SVM accuracy at each sample
for i_sample = 1:nSamples

    [p, h, stats] = signrank(svm_acc_time(:, i_sample), .5, 'tail', 'right');
    
    p_timecourse(i_sample)    = p;
    stat_timecourse(i_sample) = stats.signedrank;

end

% find clusters in the original data timecourse
clusters_orig = find_temporal_clusters(stat_timecourse, p_timecourse, p_thresh);

clusters_orig.inputs.svm_acc_time = svm_acc_time;

toc

%% analysis of repeatedly shuffled data set

% shuffle

% h = waitbar(0, 'Shuffling...');

for i_rep = 1:nReps

    % waitbar(i_rep / nReps);

    i_rep
    fix(clock)
    
    % tic
    
    clear svm_acc_time
    % get SVM accuracy for each subject and sample
    for i_sub = 1:nSubs

        labels_i = labels{i_sub};
        
        % shuffle trial labels for this subject/repetition
        nTrials  = length(labels_i);
        labels_shuffled = labels_i(randperm(nTrials));

        for i_sample = 1:nSamples

            data_i   = data{i_sub}{i_sample};

            svm_acc = libsvm_nfold(labels_shuffled, data_i, nfolds, nreps_svm, libsvm_options_train, libsvm_options_test, scale_data);

            svm_acc_time(i_sub, i_sample) = mean(svm_acc(:));

        end
    end

    % get statistical summary of SVM accuracy at each sample
    for i_sample = 1:nSamples

        [p, h, stats] = signrank(svm_acc_time(:, i_sample), .5, 'tail', 'right');

        p_timecourse_shuffle(i_sample)    = p;
        stat_timecourse_shuffle(i_sample) = stats.signedrank;

    end

    % find clusters in the shuffled data topo and get the max cluster stat
    clusters_shuffle{i_rep} = find_temporal_clusters(stat_timecourse_shuffle, p_timecourse_shuffle, p_thresh);
% %     shuffleMaxStat(i_rep) = clusters_shuffle{i_rep}.maxStatSumAbs;
    shuffleMaxStatPos(i_rep) = clusters_shuffle{i_rep}.maxStatSumPos;
    
    % toc

end

% close(h);

%% calculate p-values for each cluster in the original data set 
%  on the basis of the estimated null distribution from the shuffled data set

for i = 1:clusters_orig.nClusters
    pval = sum( shuffleMaxStatPos > clusters_orig.cluster_statSum(i) ) / nReps;
    clusters_orig.cluster_pval(i) = pval;
end


end