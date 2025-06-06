function [svm_acc, models] = libsvm_nfold_crossdecode(group, data_train_t, data_test_t, nfolds, nreps, libsvm_options_train, libsvm_options_test, scale_data)
% [svm_acc, models] = libsvm_nfold(group, data_train_t, data_test_t, nfolds, nreps, libsvm_options_train, libsvm_options_test, scale_data)
% 
% group        - label for each trials
% data_train_t - nTrials x nAttributes matrix, at training timepoint
% data_test_t  - nTrials x nAttributes matrix, at testing timepoint
% nfolds       - number of training/testing folds
% nreps        - number of repetitions of balanced training subsets in each fold
% libsvm_options_train - libsvm option string for SVM training
% libsvm_options_test  - libsvm option string for SVM testing
% scale_data   - if 1, data in each feature is rescaled to [-1, 1] for each
%                testing set. [default = 1]
% 
% svm_acc - svm accuracy for each fold x repetition
% models  - fold x repetition cell array holding the SVM models as output
%           by libsvmtrain for each fold and repetition 


if ~exist('scale_data', 'var') || isempty(scale_data)
    scale_data = 1;
end


ntrials = length(group);

% define fold indeces by interleaving consecutive trials
for i_fold = 1:nfolds
    ind_fold{i_fold} = i_fold : nfolds : ntrials;
end


% perform training and testing for each fold
for i_fold = 1:nfolds
    
    ind_test  = ind_fold{i_fold};
    ind_train = setdiff(1:ntrials, ind_test);

    % define training and testing sets
    data_train = data_train_t(ind_train, :);
    data_test  = data_test_t(ind_test, :);

    group_train = group(ind_train);
    group_test  = group(ind_test);


    % count number of trials for each group level to allow for svm
    % training that has balanced number of trials for each group level
    ng1 = sum(group_train ==  1);
    ng2 = sum(group_train == -1);

    if ng1 < ng2, ngmin = ng1; else ngmin = ng2; end


    % repeatedly take a random, balanced subset of training trials,
    % train the classifier, and test classification performance on the
    % test set
    parfor rand_rep = 1:nreps

        % get random balanced set for training data group 1
        data_train1  = data_train(group_train == 1, :);
        group_train1 = group_train(group_train == 1);

        ind_rand1 = randperm(length(group_train1));
        ind_bal1  = ind_rand1(1:ngmin);
        data_train1  = data_train1(ind_bal1, :);
        group_train1 = group_train1(ind_bal1);

        % get random balanced set for training data group 2
        data_train2  = data_train(group_train == -1, :);
        group_train2 = group_train(group_train == -1);

        ind_rand2 = randperm(length(group_train2));
        ind_bal2  = ind_rand2(1:ngmin);
        data_train2  = data_train2(ind_bal2, :);
        group_train2 = group_train2(ind_bal2);

        % concatenate
        data_train_i = [data_train1; data_train2];
        group_train_i = [group_train1; group_train2];

        % train and classify
        model = libsvmtrain(group_train_i, data_train_i, libsvm_options_train, scale_data);
        predicted_group = libsvmpredict(group_test, data_test, model, libsvm_options_test);        
        
        svm_acc(i_fold, rand_rep) = mean(predicted_group == group_test);
        models{i_fold, rand_rep} = model;
    end
end

end