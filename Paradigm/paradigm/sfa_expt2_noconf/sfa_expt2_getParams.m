function param = sfa_expt2_getParams(window, freq, isBB, isMEG)

%% stimulus

load series_selection_scale2Hz

param.f_sample_inHz    = freq;

param.nTonesInSeries   = length(series_Hz{1}(1,:)) + 1;
param.toneDur_inSecs   = .3;
param.seriesDur_inSecs = param.toneDur_inSecs * param.nTonesInSeries;

param.paBuffer_inSecs       = .3; % time between loading stim in PA buffer and playing it
param.fixationBuffer_inSecs = .5; %.3; % time between offset of last trial's FB and onset of this trial's fixation dot

param.FBoffset2stimOnset_inSecs = 1.2;
param.trialStartBuffer_inSecs = param.FBoffset2stimOnset_inSecs - param.paBuffer_inSecs;

%% display

param.BGcolor   = 0; 255; %127;
param.font      = 'Helvetica';
param.fontColor = 255; 0; %0;
param.textSizeTrial = 30;
param.textSizeBlock = 25;
param.textSize  = param.textSizeTrial;
param.sx        = 200;
param.vSpacing  = 1.4;

% fixation point
param.fixationRadius_inDegrees = .5;
if isMEG
    param.distFromScreen_inCm = 60;
else
    param.distFromScreen_inCm = 60;
end

[midW, midH] = getScreenMidpoint(window);
fixationRadius_inPixels = degrees2pixels(param.fixationRadius_inDegrees, param.distFromScreen_inCm);
param.fixationRect = [midW - fixationRadius_inPixels, ...
                      midH - fixationRadius_inPixels, ...
                      midW + fixationRadius_inPixels, ...
                      midH + fixationRadius_inPixels];

%% response

% input device
if IsOSX
    param.keyboardNumber = getKeyboardNumber;
else
    param.keyboardNumber = [];
end

param.isBB = isBB;


% time allowed for each response before triggering next trial automatically
param.respDur_inSecs = 5;


% input keys
param.exitKey      = 'ESCAPE';

if isBB
    param.readyKey = 'Return';
    param.r1_validKeys = {'0)' '1!' '2@' '3#' '4$' param.exitKey};
    param.r2_validKeys = {'5%' '6^' '7&' '8*' '9(' param.exitKey};
    
else
%     param.readyKey = [];
    param.readyKey = 'Return';
    param.r1_validKeys = {'1!' '2@' '3#' '4$' '5%' param.exitKey};
    param.r2_validKeys = {'8*' '9(' '0)' '-_' '=+' param.exitKey};    
    
end

param.r3_validKeys = param.r1_validKeys;





%% block structure

param.nTrialsPerBlock     = 30;  % should be divisible by 360
param.breakDur_inSecs     = Inf; % if Inf, break is ended by subject
param.breakWarning_inSecs = 10;  % point at which time left starts to count down
param.breakBuffer_inSecs  = 1;   % buffer b/t end of break and start of next trial


%% triggers

param.io_address = 16376;
param.io_address = hex2dec('3000');

param.io_address = 16376;
param.io_address = hex2dec('3FF8');

param.triggerDur_inSecs = .02;

param.trig_start = 255;
param.trig_stim  = [11 12 13 14 15];
param.trig_r1    = [21 22 23 24 25 29];
param.trig_r2    = [31 32 33 34 35 39];
param.trig_r3    = [41 42 43 44 45 49];

param.isMEG = isMEG;