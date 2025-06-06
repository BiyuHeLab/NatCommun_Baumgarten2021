function [predicted_label, accuracy, decision_values] = libsvmpredict(testing_label_vector, testing_instance_matrix, model, libsvm_options)
% [predicted_label, accuracy, decision_values/prob_estimates] = ...
%     libsvmpredict(testing_label_vector, testing_instance_matrix, model [, 'libsvm_options']);
% 
% This is a wrap-around function for libsvm's svmpredict function in Matlab.
%
% INPUTS
% ------
% -testing_label_vector:
%     An m by 1 vector of prediction labels. If labels of test
%     data are unknown, simply use any random values. (type must be double)
% -testing_instance_matrix:
%     An m by n matrix of m testing instances with n features.
%     It can be dense or sparse. (type must be double)
% -model:
%     The output of svmtrain.
% -libsvm_options:
%     A string of testing options in the same format as that of LIBSVM.
%     
% options:
% -b probability_estimates: whether to predict probability estimates, 0 or 1 (default 0); for one-class SVM only 0 is supported
% 
% 
% 
% The k in the -g option means the number of attributes in the input data.
% 
% option -v randomly splits the data into n parts and calculates cross
% validation accuracy/mean squared error on them.
% 
% See libsvm FAQ for the meaning of outputs.
% 
% 
% OUTPUTS
% -------
% - predicted_label:
%     a vector of predicted labels
% - accuracy:
%     a vector including accuracy (for classification), mean
%     squared error, and squared correlation coefficient (for regression).
% - decision_values:
%     a matrix containing decision values or probability
%     estimates (if '-b 1' is specified). If k is the number of classes
%     in training data, for decision values, each row includes results of 
%     predicting k(k-1)/2 binary-class SVMs. For classification, k = 1 is a
%     special case. Decision value +1 is returned for each testing instance,
%     instead of an empty vector. For probabilities, each row contains k values
%     indicating the probability that the testing instance is in each class.
%     Note that the order of classes here is the same as 'Label' field
%     in the model structure.

% default libsvm_options = []
if ~exist('libsvm_options', 'var')
    libsvm_options = [];
end

% scale the testing data in the same way training data was scaled
m = model.scaling.m;
b = model.scaling.b;

nTrials  = size(testing_instance_matrix, 1);
testing_instance_matrix_sc = testing_instance_matrix .* repmat(m, nTrials, 1) - repmat(b, nTrials, 1);   

% remove 'scaling' field so svmpredict can recognize model
model = rmfield(model, 'scaling');

% run the model on the testing data
[predicted_label, accuracy, decision_values] = svmpredict(testing_label_vector, testing_instance_matrix_sc, model, libsvm_options);