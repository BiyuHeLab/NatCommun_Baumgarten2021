function data_MEG = baseline_correct( data_MEG, baseline_inSec, fs )
% data_MEG = baseline_correct( data_MEG, baseline_inSec, fs )
%
% baseline correct a set of MEG data by subtracting from each sensor the
% average value in that sensor over the first "baseline_inSec" seconds
%
% format of data_MEG is nSensors x nSamples

% baseline_inSec = .05;
baseline_inSamples = baseline_inSec * fs;

nSensors = size(data_MEG, 1);
for i_sensor = 1:nSensors
    baseline = mean( data_MEG(i_sensor, 1 : baseline_inSamples) );
    data_MEG(i_sensor, :) = data_MEG(i_sensor, :) - baseline;
end
    

end