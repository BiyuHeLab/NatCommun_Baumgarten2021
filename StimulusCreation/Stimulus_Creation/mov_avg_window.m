function [smoothed, index] = mov_avg_window(series, nSamples)
% [smoothed, index] = mov_avg_window(series, nSamples)
%
% Compute a moving average window on input vector "series", where each
% window is nSamples long.
%
% "smoothed" is the result of performing the moving average window.
% "index" lists the indeces in "series" that have corresponding values in
%         "smoothed". Note that samples from the beginning and the end of
%         "series" will not have corresponding values in "smoothed".

offset = (nSamples-1)/2;

index = 1+offset : length(series) - offset;

for i = index
    smoothed(i-offset) = mean(series(i-offset:i+offset));
end