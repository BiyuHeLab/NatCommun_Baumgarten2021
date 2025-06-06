function predicted_label = matsvmpredict(testing_label_vector, testing_instance_matrix, model)
% [predicted_label, accuracy, decision_values/prob_estimates] = ...
%     matsvmpredict(testing_label_vector, testing_instance_matrix, model);
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
% -s svm_type : set type of SVM (default 0)
% 	0 -- C-SVC		(multi-class classification)
% 	1 -- nu-SVC		(multi-class classification)
% 	2 -- one-class SVM	
% 	3 -- epsilon-SVR	(regression)
% 	4 -- nu-SVR		(regression)
% -t kernel_type : set type of kernel function (default 2)
% 	0 -- linear: u'*v
% 	1 -- polynomial: (gamma*u'*v + coef0)^degree
% 	2 -- radial basis function: exp(-gamma*|u-v|^2)
% 	3 -- sigmoid: tanh(gamma*u'*v + coef0)
% 	4 -- precomputed kernel (kernel values in training_set_file)
% -d degree : set degree in kernel function (default 3)
% -g gamma : set gamma in kernel function (default 1/num_features)
% -r coef0 : set coef0 in kernel function (default 0)
% -c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)
% -n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)
% -p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)
% -m cachesize : set cache memory size in MB (default 100)
% -e epsilon : set tolerance of termination criterion (default 0.001)
% -h shrinking : whether to use the shrinking heuristics, 0 or 1 (default 1)
% -b probability_estimates : whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)
% -wi weight : set the parameter C of class i to weight*C, for C-SVC (default 1)
% -v n: n-fold cross validation mode
% -q : quiet mode (no outputs)
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
ntest  = size(testing_instance_matrix, 1);
testing_instance_matrix_sc = testing_instance_matrix .* repmat(m, ntest, 1) - repmat(b, ntest, 1);   

% remove 'scaling' field so svmpredict can recognize model
model = rmfield(model, 'scaling');

% run the model on the testing data
predicted_label = svmclassify(model, testing_instance_matrix_sc);
