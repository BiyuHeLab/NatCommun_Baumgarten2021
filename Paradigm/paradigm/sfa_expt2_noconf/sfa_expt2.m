clear

%% restart?
% if restarting, 
% save data in online backup to new folder first
restarting = 0;


%% startup

nPracticeTrials = 10;

% set up directory structure
expName = 'sfa_expt2';
addpath([pwd '/functions']);

valid = 0;
while ~valid
    isPractice = input('Is this a practice session? (y/n) ...   ','s');
    
    valid = any( strcmp(isPractice, {'y' 'n'}) );
    if ~valid
        fprintf('invalid input\n\n');
        
    else
        switch isPractice
            case 'y', isPractice = 1;
            case 'n', isPractice = 0;
        end
        
        if isPractice
            fprintf('This is a practice session of %d trials.\n\n', nPracticeTrials);            
        else
            fprintf('This is a full experimental session.\n\n');            
        end
    end
end

startup = startExperiment(expName, isPractice);

%% check for triggers and button box

valid = 0;
while ~valid
    isMEG = input('Is this an MEG session? (y/n) ...   ','s');
    
    valid = any( strcmp(isMEG, {'y' 'n'}) );
    if ~valid
        fprintf('invalid input\n\n');
        
    else
        switch isMEG
            case 'y', isMEG = 1;
            case 'n', isMEG = 0;
        end
        
        if isMEG
            fprintf('This is an MEG session. Trigger functionality is turned on.\n\n');            
        else
            fprintf('This is not an MEG session. Trigger functionality is shut off.\n\n');            
        end
    end
end


valid = 0;
while ~valid
    isBB = input('Are we using button boxes for input? (y/n) ...   ','s');
    
    valid = any( strcmp(isBB, {'y' 'n'}) );
    if ~valid
        fprintf('invalid input\n\n');
        
    else
        switch isBB
            case 'y', isBB = 1;
            case 'n', isBB = 0;
        end        
        
        if isBB
            fprintf('Input will be collected through the button boxes.\n\n');
        else
            fprintf('Input will be collected through the keyboard.\n\n');
        end
    end
end
WaitSecs(2);


%% launch PTB

% initialize PsychPortAudio
fprintf('Initializing PsychPortAudio...\n');
fprintf('------------------------------\n\n');

AssertOpenGL;

lowLatency = 1;
InitializePsychSound( lowLatency );

deviceID       = [];
mode           = 1;  % sound playback only
timingPriority = 2;  % take aggressive control of audio device
freq           = []; % allow PsychPortAudio to set the determined optimal value, depending on other settings 
nChannels      = 2;
bufferSize     = 128*2^3;
pahandle       = PsychPortAudio('Open', deviceID, mode, timingPriority, freq, nChannels, bufferSize);

s    = PsychPortAudio('GetStatus', pahandle)
freq = s.SampleRate;



if isMEG
    fprintf('Initializing port I/O...\n');
    fprintf('-------------------------\n\n');

    config_io;
end



fprintf('\n\nReady to go?\n');
proceed = 0;
while ~proceed
    control = input('y --> proceed\nd --> debug mode\nq --> quit\n\n...   ','s');
    
    valid = any( strcmp( control, {'y' 'd' 'q'} ) );
    if ~valid
        fprintf('invalid input\n\n');
        
    else
        switch control
            case 'y', proceed = 1;
            case 'q', PsychPortAudio('Close', pahandle); return;
            case 'd', keyboard;
        end
    end
end


% launch PTB window
[window keyboardNumber] = openScreen;


%% get experiment parameters

param = sfa_expt2_getParams(window, freq, isBB, isMEG);
param.restarting = restarting;

Screen('FillRect',  window, param.BGcolor);
Screen('TextSize',  window, param.textSize);
Screen('TextColor', window, param.fontColor);
Screen('TextFont',  window, param.font);

DrawFormattedText(window,'getting ready...','center','center');
Screen('Flip',window);



%% create stimuli for the main experiment
if ~restarting
    stim  = sfa_expt2_makeStim(param, pahandle);
    stim_ind_order = stim.ind_order;
    save([pwd '\data\online_backup\stim_ind_order.mat'], 'stim_ind_order');

else
    % load previously created stim if we're restarting
    load([pwd '\data\online_backup\stim_ind_order.mat']);
    param.stim_ind_order = stim_ind_order;
    stim  = sfa_expt2_makeStim(param, pahandle);
    
end



%% run the experiment

if param.isMEG
    DrawFormattedText(window,'Please wait a moment while\n\nthe experimenters get everything set up.','center','center');
    Screen('Flip',window);
    
    % wait indefinitely, to allow experimenter to get MEG computer
    % ready for next recording
    bkey = recordValidKeys(GetSecs, Inf, param.keyboardNumber, {param.readyKey param.exitKey});
    if strcmp(bkey, param.exitKey)
        exitNow = 1;
        return
    end            

    sx = 'center';
    sy='center';
    blockText = 'Press any key to continue.'; 
    DrawFormattedText(window, blockText, sx, sy, param.fontColor);
    Screen('Flip',window);

    WaitSecs(1);
    KbWait;

else            
    DrawFormattedText(window,'Ready to start.\n\nPress any key to continue.','center','center');
    Screen('Flip',window);

    WaitSecs(.5);
    KbWait;

end

param.isPractice = isPractice;
if isPractice
    param.nTrialsPerBlock = nPracticeTrials;
    stim.nTrials = nPracticeTrials;
end
    
[data exitNow] = sfa_expt2_runBlock(window, param, stim);

Screen('CloseAll');
PsychPortAudio('Close', pahandle);

startup.endTime    = clock;
startup.endTimeStr = datestr(startup.endTime);
stim.soundwave  = [];
save(startup.dataFile, 'startup', 'param', 'stim', 'data');

