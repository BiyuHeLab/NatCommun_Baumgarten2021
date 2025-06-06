function logL = probit_logL(x, beta, mu, y)
% logL = probit_logL(x, beta, mu, y)
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

[nTrials nX] = size(x);

logL = 0;
for i = 1:nTrials
    prob = probit_p( x(i,:), beta, mu ); % get the full set of model-predicted probabilities for this x
    p    = prob(y(i) + 1);  % get the particular probability of the outcome y for this trial
    logL = logL + log( p ); % increment the log likelihood (implicit multinomial model)
end