function [series_filtered, gain, freq] = butterworth_biyu(series, fs, f_cut, order, type, plot_gain_sq)

% [series_filtered, gain, freq] = butterworth_biyu(series, fs, f_cut, order, type, plot_gain)
% 
% Filter a time series using a Butterworth filter, following the procedures outlined
% in the supplementary materials of He et al 2010 Neuron.
% 
% Inputs
% ------
% series - The time series to be filtered. If this has an odd number of samples N+1, 
%          then only the first N samples will be used to create the filtered series.
% 
% fs     - Sampling frequency
% 
% f_cut  - The cutoff frequency for filtering. Should contain only one value for lowpass
%          and highpass filters. For bandpass filters, should be 1x2 vector [f1, f2], f1 < f2.
% 
% order  - Order of the Butterworth filter. Higher values have a steeper rolloff around 
%          the cutoff frequency.
%          
% type   - Type of filter. Can be 'lowpass', 'highpass', or 'bandpass'. Bandpass filters
%          for frequency range [f1, f2] are constructed by multiplying a lowpass filter
%          for f2 with a highpass filter for f1.   [default = 'lowpass']
%              
% plot_gain_sq - Make a plot of the squared gain of the Butterworth filter as a function 
%                of frequency. Set to 1 to make the plot, 0 for no plot.  [default = 0]
% 
% Outputs
% -------
% series_filtered - The filtered time series.
% gain            - The gain of the Butterworth filter as a function of frequency.
% freq            - The frequencies for the gain values of the Butterworth filter.


%% handle inputs

% set default values for variables
if ~exist('type','var') || isempty(type)
    type = 'lowpass';
end

if ~exist('plot_gain_sq','var') || isempty(plot_gain_sq)
    plot_gain_sq = 0;
end


% check inputs
switch type
    case 'lowpass'
        if length(f_cut) > 1
            error('for low pass filter, use only one cutoff frequency')
        end
        
        
    case 'highpass'
        if length(f_cut) > 1
            error('for high pass filter, use only one cutoff frequency')
        end
        
        
    case 'bandpass'
        if length(f_cut) ~= 2
            error('for bandpass filter, specify two cutoff frequencies')
        elseif f_cut(1) >= f_cut(2)
            error('for bandpass filter, first cutoff frequency must be lower than the second')
        end
end


% ensure series has even number of samples
if rem(length(series), 2) == 1
    series = series(1:end-1);
end


%% set up variables

nSamples = length(series);

% number of seconds in the series
period = nSamples / fs;

% increments in Hz on frequency domain plot
hzpbin = 1 / period;

freq = 0 : hzpbin : (nSamples-1) * hzpbin;


%% compute FFT

fft_series = fft( detrend(series) );


%% compute the butterworth filter weights

switch type
    case 'lowpass'
        r    = ( freq ./ f_cut ).^(2*order);
        gain = sqrt( 1 ./ (1 + r) );
        
    case 'highpass'
        r    = ( freq ./ f_cut ).^(2*order);
        gain = sqrt( r ./ (1 + r) );
       
    case 'bandpass'
        % for bandpass in frequency band [f1, f2],
        % we multiply the lowpass filter for f2 and the highpass filter for f1
        r_low    = ( freq ./ f_cut(2) ).^(2*order);
        gain_low = sqrt( 1 ./ (1 + r_low) );        

        r_high    = ( freq ./ f_cut(1) ).^(2*order);
        gain_high = sqrt( r_high ./ (1 + r_high) );        
        
        gain = gain_low .* gain_high;

end

% mirror the gain function above and below nyquist frequency
g1 = gain(1:nSamples/2+1);
g2 = gain(nSamples/2 : -1 : 2);  
gain = [g1 g2];


%% plot the gain of the Butterworth filter

if plot_gain_sq
    switch type
        case 'bandpass'

            g1 = gain_low(1:nSamples/2+1);
            g2 = gain_low(nSamples/2 : -1 : 2);  
            gain_low = [g1 g2];

            g1 = gain_high(1:nSamples/2+1);
            g2 = gain_high(nSamples/2 : -1 : 2);  
            gain_high = [g1 g2];        

            figure; grid on; hold on;
            plot(freq, gain_low.^2, 'r-');
            plot(freq, gain_high.^2, 'g-');
            plot(freq, gain.^2, 'k-');
            xlim([0 fs/2]);
            xlabel('frequency in Hz')
            ylabel('power (a.u.)')

        otherwise
            figure; grid on; hold on;
            plot(freq, gain.^2, 'b-');
            xlim([0 fs/2]);            
            xlabel('frequency in Hz')
            ylabel('power (a.u.)')
            
    end
end

%% apply the filter

fft_filtered    = fft_series .* gain;
series_filtered = ifft(fft_filtered);