function [freq, power, phase, amplitude] = fft_window_dataset(data, trials, sensors, fs, windowDur_inSecs, windowOverlap_inSecs, hanning, freq_cutoff)

% [freq, power, phase, amplitude] = fft_window_dataset(data, trials, sensors, fs, windowDur_inSecs, windowOverlap_inSecs, hanning, freq_cutoff)
% 
% helper function which computes FFT using fft_window for a specified set
% of trials and sensors in the Fieldtrip "data" struct
%
% INPUTS
%
% data    - the EEG/MEG data to be analyzed, as formatted by fieldtrip
%           specifically, should have a cell array field data.trial of length
%           nTrials, where each element is a matrix with nSensors rows and
%           nSamples columns
% trials  - the indeces to data.trial which should be extracted for FFT
%           averaging
% sensors - the indeces the sensors (rows) of each data.trial{} matrix
%           which should be extracted for FFT averaging
%
% All remaining inputs are as described in the help for the function fft_window
%
%
%
% OUTPUTS
%
% freq - contains the frequency for each element of the FFT analysis
% power, phase, amplitude - contain FFT outputs at each frequency for each
%        sensor, averaged across trials and FFT windows. matrices with
%        length(sensors) rows and length(freq) columns


power_forEachSensor = [];
phase_forEachSensor = [];
amp_forEachSensor   = [];

for i_sensor = 1:length(sensors)
    
    s = sensors(i_sensor);

    power_over_trials = [];
    phase_over_trials = [];
    amp_over_trials   = [];    
    
    for i_trial = 1:length(trials)

        t = trials(i_trial);
        
        % do FFT analysis for current trial / sensor
        series = data.trial{t}(s, :);
        [freq, power, phase, amplitude] = fft_window(series, fs, windowDur_inSecs,  windowOverlap_inSecs, hanning);

            
        % average across windows, if applicable
        if size(power,1) > 1
            power     = mean(power);
            phase     = mean(phase);
            amplitude = mean(amplitude);
        end

        
        % exclude data for high frequencies
        f     = freq <= freq_cutoff;

        freq      = freq(f);
        power     = power(f);
        phase     = phase(f);
        amplitude = amplitude(f);

        
        % store results for analysis of current trial
        power_over_trials(i_trial, :) = power;
        phase_over_trials(i_trial, :) = phase;
        amp_over_trials(i_trial, :)   = amplitude;
        
    end
    
    
    % average across trials, if applicable
    if length(trials) > 1
        power_over_trials = mean(power_over_trials);
        phase_over_trials = mean(phase_over_trials);
        amp_over_trials   = mean(amp_over_trials);
    end
    
    power_forEachSensor(i_sensor, :) = power_over_trials;
    phase_forEachSensor(i_sensor, :) = phase_over_trials;
    amp_forEachSensor(i_sensor, :)   = amp_over_trials;
    
end


% simplify output variable names
power     = power_forEachSensor;
phase     = phase_forEachSensor;
amplitude = amp_forEachSensor;