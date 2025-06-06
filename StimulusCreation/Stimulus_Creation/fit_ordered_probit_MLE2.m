function [mu sd logL] = fit_ordered_probit_MLE2(x, y, m)

mu_guess = linspace(min(x), max(x), m+2);
mu_guess = mu_guess(2:end-1);

sd_guess = 1;

guess = [mu_guess sd_guess];


[N V] = size(x);

A = []; B = [];
for i = 1:m-1
    nz1 = i-1;
    nz2 = length(guess) - (2+nz1);
    A = [A; zeros(1,nz1) 1 -1 zeros(1,nz2)];
    B = [B; -1e-5];
end

LB = [zeros(size(mu_guess)) 1e-10];

    
save('fit.mat','x','y','m');
% global x y m

% [param nlogL] = fminunc(@probit_nlogL, guess);
[param nlogL] = fmincon(@probit_nlogL, guess, A, B, [], [], LB);


mu   = param(1:m);
sd   = param(m+1);
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
sd   = param(m+1);

% if ~isreal(mu) || ~isreal(sd)
%     keyboard
% end

logL = 0;
for i = 1:N
    prob = probit_p( x(i,:), mu, sd ); % get the full set of model-predicted probabilities for this x
    p    = prob(y(i) + 1);  % get the particular probability of the outcome y for this trial
    logL = logL + log( p ); % increment the log likelihood (implicit multinomial model)
end

if ~isreal(logL)
    keyboard
end

if isnan(logL)
    nlogL = Inf;
else
    nlogL = -logL;
end

end


function prob = probit_p(x, mu, sd)
% prob = probit_p(x, beta, mu)
%
% Given an independent variable x, weighting factor beta, and thresholds mu, 
% calculate the probability of ordinal output variable 
%
% y = j, for j = 0 to m.
%
% where y* = x*beta + e, and y = j iff y* is located between the
% appropriate thresholds in mu.
%
% Inputs
% ------
% x    - 1 x V vector of independent variables
% beta - V x 1 vector of weights
% mu   - 1 x m vector of thresholds
%
% Output
% ------
% prob - 1 x m+1 vector of probabilities. 
%        For each entry i, prob(i) = p( y == i-1 | x, beta, mu )


mu = [-Inf mu Inf];
m  = length(mu) - 2;

prob = normcdf( mu(2:end) - x, 0, sd ) - normcdf( mu(1:end-1) - x, 0, sd );

end