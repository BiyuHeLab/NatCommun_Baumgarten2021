function model = matsvmtrain(training_label_vector, training_instance_matrix, scale_data)
% model = matsvmtrain(training_label_vector, training_instance_matrix [, scale_data]);
% 
% This is a wrap-around function for libsvm's svmtrain function in Matlab.
%
% INPUTS
% ------
% -training_label_vector:
%     An m by 1 vector of training labels (type must be double).
% -training_instance_matrix:
%     An m by n matrix of m training instances with n features.
%     It can be dense or sparse (type must be double).
% -libsvm_options:
%     A string of training options in the same format as that of LIBSVM.
% - scale_data
%     if 1, data in each feature is rescaled to [-1, 1]. [default = 1]
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
% The 'svmtrain' function returns a model which can be used for future
% prediction.  It is a structure and is organized as [Parameters, nr_class,
% totalSV, rho, Label, ProbA, ProbB, nSV, sv_coef, SVs]:
% 
%         -Parameters: parameters
%         -nr_class: number of classes; = 2 for regression/one-class svm
%         -totalSV: total #SV
%         -rho: -b of the decision function(s) wx+b
%         -Label: label of each class; empty for regression/one-class SVM
%         -sv_indices: values in [1,...,num_traning_data] to indicate SVs in the training set
%         -ProbA: pairwise probability information; empty if -b 0 or in one-class SVM
%         -ProbB: pairwise probability information; empty if -b 0 or in one-class SVM
%         -nSV: number of SVs for each class; empty for regression/one-class SVM
%         -sv_coef: coefficients for SVs in decision functions
%         -SVs: support vectors
% 
% If you do not use the option '-b 1', ProbA and ProbB are empty
% matrices. If the '-v' option is specified, cross validation is
% conducted and the returned model is just a scalar: cross-validation
% accuracy for classification and mean-squared error for regression.
% 
% More details about this model can be found in LIBSVM FAQ
% (http://www.csie.ntu.edu.tw/~cjlin/libsvm/faq.html) and LIBSVM
% implementation document
% (http://www.csie.ntu.edu.tw/~cjlin/papers/libsvm.pdf).


if ~exist('scale_data', 'var') || isempty(scale_data)
    scale_data = 1;
end


%% scale the data

nTrials   = size(training_instance_matrix, 1);
nFeatures = size(training_instance_matrix, 2);

if scale_data
    % scale the training data such that each attribute ranges from [-1, 1]
    drange    = range(training_instance_matrix);
    m = 2 ./ drange;
    b = (m .* max(training_instance_matrix)) - 1;

else
    m = ones(1, nFeatures);
    b = zeros(1, nFeatures);
    
end

training_instance_matrix_sc = training_instance_matrix .* repmat(m, nTrials, 1) - repmat(b, nTrials, 1);



%% run the svm training
w = cd('/usr/local/MATLAB/R2015a/toolbox/stats/stats');
opt = optimset;
opt.MaxIter = 15000 * 1000;
model = svmtrain(training_instance_matrix_sc, training_label_vector, 'method', 'SMO', 'options', opt);
cd(w); 

        
% save the data scaling parameters for use with the testing data
model.scaling.m = m;
model.scaling.b = b;