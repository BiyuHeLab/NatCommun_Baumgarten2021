function [data exitNow] = sfa_expt2_runTrial(window, param, stim, trialNum, data)

data.timing.trialStart(trialNum) = GetSecs;

exitNow = 0;

data.trialNum(trialNum) = trialNum;

sx = 'center';
sy = 'center';
vSpacing = 1.4;


%% show fixation dot

% give some buffer time between the offset of last trial's FB and the onset
% of this trial's fixation
WaitSecs(param.fixationBuffer_inSecs);

Screen('FillOval', window, param.fontColor, param.fixationRect);
data.timing.fixationOnset(trialNum) = Screen('Flip', window);

% % % image = Screen ('GetImage', window);
% % % file = 'sfa2_fixation.png';
% % % imwrite(image,file,'png');


%% buffer time so that auditory stim onset is 1 sec after FB offset of last trial

if trialNum == 1
    WaitSecs(param.trialStartBuffer_inSecs - param.fixationBuffer_inSecs);
else
    waitUntil(data.timing.FBoffset(trialNum-1), param.trialStartBuffer_inSecs);   
end


%% play the stimulus

t0 = GetSecs + param.paBuffer_inSecs;

% trial_text = ['trial ' num2str(trialNum) ' of 360'];
% DrawFormattedText(window, trial_text, sx, sy, param.fontColor, [], [], [], vSpacing);
% Screen('Flip', window);

% load buffer
w  = stim.soundwave{trialNum};
w2 = [w; w];

PsychPortAudio('FillBuffer', stim.pahandle, w2);


% initialize playback
data.timing.stimOnset(trialNum)     = PsychPortAudio('Start', stim.pahandle, [], t0);
data.timing.stimOnset_req(trialNum) = t0;

% send the trigger
if param.isMEG
    trigger = param.trig_stim( stim.betaID(trialNum) );
    data.triggers.stim_trigger(trialNum) = trigger;
    data.timing.stim_trigger(trialNum)   = send_trigger(trigger, param);
end

% check pahandle status until it's done
WaitSecs(param.paBuffer_inSecs + param.seriesDur_inSecs - 2);

s = PsychPortAudio('GetStatus', stim.pahandle);
while s.Active == 1
    s = PsychPortAudio('GetStatus', stim.pahandle);
end

% % remove fixation
% Screen('Flip', window);

% estimate durtaion of series playback with two methods
data.timing.dur_loop(trialNum) = GetSecs - t0;
data.timing.dur_est(trialNum)  = s.EstimatedStopTime - t0;
data.pahandle_status(trialNum) = s;


%% get behavioral responses

%%%%%%%%%%% start probability rating %%%%%%%%%%%

% prompt for final tone probability rating
% % % r1_text = 'Final tone likelihood?\n\n1 = very unlikely,  2 = somewhat unlikely,  3 = somewhat likely,  4 = very likely\n(use LEFT HAND)';
r1_text = 'Final tone likelihood?\n\n0 = very unlikely   4 = very likely\n(use LEFT HAND)';

DrawFormattedText(window, r1_text, sx, sy, param.fontColor, [], [], [], vSpacing);

WaitSecs(.4);
data.timing.r1_prompt(trialNum) = Screen('Flip',window);

% % % image = Screen ('GetImage', window);
% % % file = 'sfa2_q1.png';
% % % imwrite(image,file,'png');


% collect tone probability rating
[r1_key data.r1_RT(trialNum) data.timing.r1_time(trialNum)] = ... 
    recordValidKeys(data.timing.r1_prompt(trialNum), param.respDur_inSecs, param.keyboardNumber, param.r1_validKeys);

Screen('Flip',window);

% assess response
id = 6;
rk = param.r1_validKeys;
switch r1_key
    case rk{1},         data.resp_prob(trialNum) =  1;  id = 1;
    case rk{2},         data.resp_prob(trialNum) =  2;  id = 2;
    case rk{3},         data.resp_prob(trialNum) =  3;  id = 3;
    case rk{4},         data.resp_prob(trialNum) =  4;  id = 4;
    case rk{5},         data.resp_prob(trialNum) =  5;  id = 5;
    case 'noanswer',    data.resp_prob(trialNum) = -1;
    case 'invalid',     data.resp_prob(trialNum) = -2;
    case 'cell',        data.resp_prob(trialNum) = -3;
    case param.exitKey, data.resp_prob(trialNum) = -4; exitNow = 1;
    otherwise,          data.resp_prob(trialNum) = -5;
end

% send trigger
if param.isMEG
    trigger = param.trig_r1(id);
    data.timing.r1_trigger(trialNum)   = send_trigger(trigger, param);
    data.triggers.r1_trigger(trialNum) = trigger;
end

% give 'too slow' warning if necessary
if strcmp(r1_key,'noanswer')
    data.timing.r1_tooSlow(trialNum) = 1;
    DrawFormattedText(window, 'Too slow!', sx, sy, param.fontColor)
    Screen('Flip', window);
    WaitSecs(2);
end

Screen('Flip',window);


%%%%%%%%%%% start beta estimation %%%%%%%%%%%

% prompt for beta estimation
% % % r2_text = 'Trend strength?\n\n7 = no trend,  8 = weak trend,  9 = strong trend\n(use RIGHT HAND)';
% r2_text = 'beta?\n\n0   --> press 1\n0.5 --> press 2\n1   --> press 3\n1.5 --> press 4\n2   --> press 5';
r2_text = 'Trend strength?\n\n0 = no trend   4 = strong trend\n(use RIGHT HAND)';

DrawFormattedText(window, r2_text, sx, sy, param.fontColor, [], [], [], vSpacing);

WaitSecs(.4);
data.timing.r2_prompt(trialNum) = Screen('Flip',window);

% % % image = Screen ('GetImage', window);
% % % file = 'sfa2_q2.png';
% % % imwrite(image,file,'png');

% collect beta estimation
[r2_key data.r2_RT(trialNum) data.timing.r2_time(trialNum)] = ... 
    recordValidKeys(data.timing.r2_prompt(trialNum), param.respDur_inSecs, param.keyboardNumber, param.r2_validKeys);

Screen('Flip',window);

% assess response
id = 6;
rk = param.r2_validKeys;
switch r2_key
    case rk{1},         data.resp_beta(trialNum) =  0;     id = 1;   tr_str_resp = 0;
    case rk{2},         data.resp_beta(trialNum) =  .5;    id = 2;   tr_str_resp = 1;
    case rk{3},         data.resp_beta(trialNum) =  1.01;  id = 3;   tr_str_resp = 2;
    case rk{4},         data.resp_beta(trialNum) =  1.5;   id = 4;   tr_str_resp = 3;
    case rk{5},         data.resp_beta(trialNum) =  2;     id = 5;   tr_str_resp = 4;
    case 'noanswer',    data.resp_beta(trialNum) = -1;
    case 'invalid',     data.resp_beta(trialNum) = -2;
    case 'cell',        data.resp_beta(trialNum) = -3;
    case param.exitKey, data.resp_beta(trialNum) = -4; exitNow = 1;
    otherwise,          data.resp_beta(trialNum) = -5;
end

% send trigger
if param.isMEG
    trigger = param.trig_r2(id);
    data.timing.r2_trigger(trialNum)   = send_trigger(trigger, param);
    data.triggers.r2_trigger(trialNum) = trigger;
end

% give 'too slow' warning if necessary
if strcmp(r2_key,'noanswer')
    data.timing.r2_tooSlow(trialNum) = 1;    
    DrawFormattedText(window, 'Too slow!', sx, sy, param.fontColor)
    Screen('Flip', window);
    WaitSecs(2);
end

Screen('Flip',window);



%%%%%%%%%%% start beta confidence rating %%%%%%%%%%%

% % % if data.timing.r2_tooSlow(trialNum) == 0
% % % 
% % % % prompt for confidence
% % % % r3_text = 'Confidence?\n\n1 = lowest   5 = highest';
% % % r3_text = 'Confidence about trend strength?\n\n0 = guess   4 = high conf\n(use LEFT HAND)';
% % % DrawFormattedText(window, r3_text, sx, sy, param.fontColor, [], [], [], vSpacing);
% % % 
% % % WaitSecs(.4);
% % % data.timing.r3_prompt(trialNum) = Screen('Flip',window);
% % % 
% % % 
% % % % collect beta estimation
% % % [r3_key data.r3_RT(trialNum) data.timing.r3_time(trialNum)] = ... 
% % %     recordValidKeys(data.timing.r3_prompt(trialNum), param.respDur_inSecs, param.keyboardNumber, param.r3_validKeys);
% % % 
% % % Screen('Flip',window);
% % % 
% % % % assess response
% % % id = 6;
% % % rk = param.r3_validKeys;
% % % switch r3_key
% % %     case rk{1},         data.conf_beta(trialNum) =  1;  id = 1;
% % %     case rk{2},         data.conf_beta(trialNum) =  2;  id = 2;
% % %     case rk{3},         data.conf_beta(trialNum) =  3;  id = 3;
% % %     case rk{4},         data.conf_beta(trialNum) =  4;  id = 4;
% % %     case rk{5},         data.conf_beta(trialNum) =  5;  id = 5;
% % %     case 'noanswer',    data.conf_beta(trialNum) = -1;
% % %     case 'invalid',     data.conf_beta(trialNum) = -2;
% % %     case 'cell',        data.conf_beta(trialNum) = -3;
% % %     case param.exitKey, data.conf_beta(trialNum) = -4; exitNow = 1;
% % %     otherwise,          data.conf_beta(trialNum) = -5;
% % % end
% % % 
% % % % send trigger
% % % if param.isMEG
% % %     trigger = param.trig_r3(id);
% % %     data.timing.r3_trigger(trialNum)   = send_trigger(trigger, param);
% % %     data.triggers.r3_trigger(trialNum) = trigger;
% % % end
% % % 
% % % % give 'too slow' warning if necessary
% % % if strcmp(r3_key,'noanswer')
% % %     data.timing.r3_tooSlow(trialNum) = 1;    
% % %     DrawFormattedText(window, 'Too slow!', sx, sy, param.fontColor)
% % %     Screen('Flip', window);
% % %     WaitSecs(2);
% % % end
% % % 
% % % Screen('Flip',window);
% % % 
% % % end



%% sort out response 


% sort out correctness
if data.resp_beta(trialNum) == stim.beta(trialNum)
    data.correct_beta(trialNum) = 1;

% N/A if response < 1
elseif data.resp_beta(trialNum) < 0
    data.correct_beta(trialNum) = -1;

% otherwise, just incorrect
else
    data.correct_beta(trialNum) = 0;
end

data.diff_beta(trialNum) = data.resp_beta(trialNum) - stim.beta(trialNum);

if exitNow, return; end


%% show feedback
tr_str_stim = stim.betaID(trialNum) - 1;

switch data.correct_beta(trialNum)
    case 1,  FB_text = ['Correct!\n\ntrend strength was ' num2str(tr_str_stim)];
    case 0,
        if abs(data.diff_beta(trialNum)) < .6
            fbt = 'Close!';
        else
            fbt = 'Incorrect!';
        end
        FB_text = [fbt '\n\ntrend strength was ' num2str(tr_str_stim) '\nyou responded trend strength = ' num2str(tr_str_resp)];
    case -1, FB_text = ['No response for trend strength recorded.\n\ntrend strength was ' num2str(tr_str_stim)];
end
DrawFormattedText(window, FB_text, sx, sy, param.fontColor, [], [], [], vSpacing);
Screen('Flip', window);
WaitSecs(2);

data.timing.FBoffset(trialNum) = Screen('Flip', window);

data.timing.trialEnd(trialNum) = GetSecs;
data.timing.trialDur(trialNum) = data.timing.trialEnd(trialNum) - data.timing.trialStart(trialNum);

end