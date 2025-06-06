function [reg, x_sample, y_sample] = reg_log_linspace(y, x, logx_min, logx_max, logx_grain)
% [reg, x_sample, y_sample] = reg_log_linspace(y, x, logx_min, logx_max, logx_grain)
% 
% Regress log10(y) onto log10(x) such that the log10(x) values are 
% approximately evenly spaced on the log10 scale.
%
% Specifically, for each value v in logx_min : logx_grain : logx_max, a
% sample from y and x is selected from the index where log10(x) is closest 
% to v 
%
% Inputs
% ------
% y, x       - the data for regression
% logx_min   - the minimum value for log10(x) used in the regression
% logx_max   - the maximum value for log10(x) used in the regression
% logx_grain - the grain utilized for the log10(x) scale
%
% Outputs
% -------
% reg      - the output from regstats(log10(y_sample), log10(x_sample))
%            see "help regstats" for more info
% x_sample - samples from x used for the regression
% y_sample - samples from y used for the regression

logx_linspace = logx_min : logx_grain : logx_max;
for i_x = 1:length(logx_linspace)
    [m, ind] = min( abs( log10(x) - logx_linspace(i_x) ) );
    x_sample(i_x) = x(ind);
    y_sample(i_x) = y(ind);
end

reg = regstats(log10(y_sample), log10(x_sample));