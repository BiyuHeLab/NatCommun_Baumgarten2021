function [beta mu logL] = fit_ordered_probit_MLE(x, y, m)

% keyboard

beta_guess = inv(x'*x)*x'*y; % OLS solution
ystar      = x * beta_guess;
mu_guess   = linspace(min(ystar), max(ystar), m+2);
mu_guess   = mu_guess(2:end-1);

guess = [mu_guess beta_guess'];


[N V] = size(x);

A = []; B = [];
for i = 1:m-1
    nz1 = i-1;
    nz2 = length(guess) - (2+nz1);
    A = [A; zeros(1,nz1) 1 -1 zeros(1,nz2)];
    B = [B; 0];
end

LB = [zeros(size(mu_guess)) -Inf*ones(size(beta_guess'))];

    
save('fit.mat','x','y','m');
% global x y m

% [param nlogL] = fminunc(@probit_nlogL, guess);
[param nlogL] = fmincon(@probit_nlogL, guess, A, B, [], [], LB);


mu   = param(1:m);
beta = param(m+1:end)';
logL = -nlogL;

% clear x y m
delete('fit.mat');

end

function nlogL = probit_nlogL(param)
% logL = probit_logL(param)
%
% Given a series of independent variables x with weighting factors beta,
% thresholds mu, and ordinal outcomes y, find the log likelihood of the
% data set in y.
%
% Inputs
% ------
% x    - N x V vector of independent variables
% beta - V x 1 vector of weights
% mu   - 1 x m vector of thresholds
% y    - 1 x N vector of ordinal outcomes
%
% Output
% ------
% logL - log likelihood of the data y, given input x and parameters beta, mu

load('fit.mat');
% global x y m

[N V] = size(x);

mu   = param(1:m);
beta = param(m+1:end)';
% beta = 1;

logL = 0;
for i = 1:N
    prob = probit_p( x(i,:), beta, mu ); % get the full set of model-predicted probabilities for this x
    p    = prob(y(i) + 1);  % get the particular probability of the outcome y for this trial
    logL = logL + log( p ); % increment the log likelihood (implicit multinomial model)
end

nlogL = -logL;

end