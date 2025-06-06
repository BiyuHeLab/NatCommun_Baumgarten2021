function [freq, power, phase, amplitude, time] = fft_window(series, fs, windowDur_inSecs,  windowOverlap_inSecs, hanning, nWindows)
% [freq, power, phase, amplitude] = fft_window(series, fs, windowDur_inSecs,  windowOverlap_inSecs, hanning, nWindows)
%
% Do FFT analysis on a time series using sliding windows.
%
% inputs
% ------
% series - the data to be analyzed. should be a vector
% fs     - sampling frequency in Hz
%
% windowDur_inSecs - duration of FFT window, in seconds
%            [ default = length(series) / fs ]
%
% windowOverlap_inSecs - extent of overlap between consecutive windows.
%    for instance, if windowDur_inSecs = 300 and windowOverlap_inSecs = 100,
%    then the first 3 FFT windows are [0, 300], [200, 500], [400, 700]
%            [ default = 0 ]
%
% hanning  - if 1, time series in the window is multiplied by a Hanning filter; 
%            if 0, no filter applied
%            [ default = 0 ]
%
% nWindows - manually enter number of windows for which to compute FFT
%            [ default = max # possible for given series ]
%
% outputs
% -------
% freq      - frequency values for FFT output
% power     - power as a function of frequency for each time window. 
%             size  = [number of windows, number of samples in series]
% phase     - likewise for phase
% amplitude - likewise for amplitude


%% handle inputs

if ~exist('windowDur_inSecs','var')
    windowDur_inSecs = length(series) / fs;
    disp(['Window duration unspecified; automatically set to length(series)/fs = ' num2str(windowDur_inSecs) ' sec']);
end

if ~exist('windowOverlap_inSecs','var')
    windowOverlap_inSecs = 0;
    disp(['Window overlap unspecified; automatically set to ' num2str(windowOverlap_inSecs) ' sec']);
end

if windowOverlap_inSecs >= windowDur_inSecs
    error('Window overlap must be less than window size');
end

if ~exist('hanning','var')
    hanning = 0;
    disp(['Hanning filter usage unspecified; automatically set to ' num2str(hanning)]);
end


%% FFT preliminaries

m = length(series);
t = (1:m) / fs;

n_window  = fs * windowDur_inSecs;
n_overlap = fs * windowOverlap_inSecs;
n_fft     = n_window; % 2^nextpow2(n);

fN   = fs/2;                       % Nyquist frequency
freq = linspace(0, fN, n_fft/2+1); % FFT freq spans from 0 to fN


%% determine the number of windows that fit into the series, 
%  and the time values on which they are centered


if ~exist('nWindows', 'var') || isempty(nWindows)
% compute maximum # of possible windows

    nWindows = 0;
    tf = 0;
    time = [];
    while tf < t(end)
        if tf == 0
            tf = tf + windowDur_inSecs;
        else
            tf = tf + windowDur_inSecs - windowOverlap_inSecs;
        end

        if tf <= t(end)
            nWindows = nWindows + 1;
            time(end+1, 1) = tf - (windowDur_inSecs/2);
        end
    end

    
else
% use # of windows specified

    tf = 0;
    time = [];
    for i_w = 1:nWindows
        if tf == 0
            tf = tf + windowDur_inSecs;
        else
            tf = tf + windowDur_inSecs - windowOverlap_inSecs;
        end
        
        time(i_w, 1) = tf - (windowDur_inSecs/2);
    end
    
end

%% calculate FFT for each window

n_fft_pos = n_fft/2 + 1;

amplitude = zeros(nWindows, n_fft_pos);
power     = zeros(nWindows, n_fft_pos);
phase     = zeros(nWindows, n_fft_pos);

for i_w = 1 : nWindows

    % determine window indeces
    if i_w == 1
        ind_min = 1;
        ind_max = ind_min + n_window - 1;
    else
        ind_min = ind_max - n_overlap + 1;
        ind_max = ind_min + n_window - 1;
    end
    
    % extract and preprocess series for this window
    s = series( ind_min : ind_max );
%     s = detrend(s);
%     s = s - mean(s);
    if hanning
        s = s .* hann(n_window)'; 
    end
    
    % FFT analysis
    y = fft(s, n_fft);
    
    % select only frequencies ranging from 0 to fN
    y = y(1 : n_fft_pos);

    % amp and power each get multiplied by 2 to reflect the "folding over"
    % of the symmetric FFT data for freq > fN
    amplitude(i_w,:) = 2 * (abs(y) / n_fft);
    power(i_w,:)     = 2 * (abs(y) / n_fft).^2;
    phase(i_w,:)     = angle(y);
    
end