function [firstEmptySample data] = find_empty_MEG_data(data)
% [firstEmptySample data] = find_empty_MEG_data(data)
% 
% Given a dataset "data" in fieldtrip format, find the sample that commences 
% a sequence of zeros at the end of the continuous block of recording stored 
% in data.trial{1}. Such a sequence of zeros may occur if the recording is
% manually ended before the pre-determined recording time is up. The index of
% the first empty sample is returned in the output "firstEmptySample".
% 
% The second output "data" contains a copy of the input data, but with all 
% samples including and after "firstEmptySample" removed from data.time{1} 
% and data.trial{1}.

[nChannels nSamples] = size(data.trial{1});

if nChannels > 1
    emptySample = all(data.trial{1} == 0);
else
    emptySample = data.trial{1} == 0;
end

firstEmptySample = nSamples + 1;
isEmpty = 1;
while isEmpty
    if emptySample( firstEmptySample-1 )
        firstEmptySample = firstEmptySample - 1;
    else
        isEmpty = 0;
    end
end

if firstEmptySample <= nSamples
    time2{1} = data.time{1}(1:firstEmptySample-1);
    trial2{1} = data.trial{1}(:, 1:firstEmptySample-1);

    data.time = time2;
    data.trial = trial2;
end