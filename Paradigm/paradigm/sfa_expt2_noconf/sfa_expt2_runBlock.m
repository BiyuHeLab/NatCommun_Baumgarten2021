function [data exitNow] = sfa_expt2_runBlock(window, param, stim)


%% initialize data struct

if ~param.restarting

% store behavioral data here
dataFields = {'trialNum' ...
    'resp_prob' 'resp_beta' 'conf_beta' 'correct_beta' 'diff_beta' ...    
    'r1_RT' 'r2_RT' 'r3_RT'}; 

dataTriggerFields = {'stim_trigger' 'r1_trigger' 'r2_trigger' 'r3_trigger'};

dataTimingFields = { ...
    'trialStart' 'trialEnd' 'trialDur' ...
    'stimOnset' 'stimOnset_req' 'stim_trigger' 'fixationOnset' ...
    'dur_loop' 'dur_est' ...
    'r1_prompt' 'r1_time' 'r1_tooSlow' 'r1_trigger' ...
    'r2_prompt' 'r2_time' 'r2_tooSlow' 'r2_trigger' ...
    'r3_prompt' 'r3_time' 'r3_tooSlow' 'r3_trigger' ...
    'FBoffset'};

nil = zeros(1,stim.nTrials);
for j = 1:length(dataFields)
    eval(['data.' dataFields{j} ' = nil;']);
end

data.triggers.block_trigger = [];
for j = 1:length(dataTriggerFields)
    eval(['data.triggers.' dataTriggerFields{j} ' = nil;']);
end

data.timing.block_trigger = [];
for j = 1:length(dataTimingFields)
    eval(['data.timing.' dataTimingFields{j} ' = nil;']);
end

s = PsychPortAudio('GetStatus', stim.pahandle);
for j = 1:stim.nTrials
    data.pahandle_status(j) = s;
end


data.stim  = stim;
data.stim.soundwave = [];
data.param = param;

else

%% but if we're restarting...

% load previously created data struct if we're restarting
load([pwd '\data\online_backup\data_sofar.mat']);


end

%% start the block

sx='center';
sy='center';

if ~param.restarting   
    trialNum = 0;
    restart_flag = 0;
else
    restart_flag = 1;
end

while trialNum < stim.nTrials
    
%     blockNum = blockNum + 1;
    blockNum = floor(trialNum/param.nTrialsPerBlock) + 1;
       
%     blockText = 'Ready for next block.\n\nPress any key to continue.';
    blockText = ['Ready for next block.\n\n' ...
                 'Remember:\n\n' ...
                 'For final tone likelihood question\n' ...
                 '- use LEFT HAND\n' ...
                 '- pinky = 0 (very unlikely), thumb = 4 (very likely)\n\n' ...
                 'For trend strength question\n' ...
                 '- use RIGHT HAND\n' ...
                 '- thumb = 0 (no trend), pinky = 4 (strong trend)\n\n' ...
                 'Press any key to continue.'];

    wrect = Screen('Rect',window);
    sx = round(wrect(3)*.1);
    Screen('TextSize',  window, param.textSizeBlock);
    DrawFormattedText(window, blockText, sx, sy, param.fontColor,[],[],[],param.vSpacing);
    Screen('Flip',window);
    Screen('TextSize',  window, param.textSizeTrial);
    
    WaitSecs(.5);
    bkey = recordValidKeys(GetSecs, Inf, param.keyboardNumber, []);
    if strcmp(bkey, param.exitKey)
        exitNow = 1;
        return
    end
    
    Screen('Flip',window);
    
    data.timing.blockStart{blockNum} = clock;
    t0 = GetSecs;
    
    % send block start trigger
    if param.isMEG
        data.timing.block_trigger(blockNum)   = send_trigger(param.trig_start, param);
        data.triggers.block_trigger(blockNum) = param.trig_start;    
    end
    
    % start at the proper point in the block, if we're restarting from a
    % previous point
    if restart_flag
        t_start = rem(trialNum, param.nTrialsPerBlock) + 1;        
        restart_flag = 0;
    else
        t_start = 1;
    end
        
    for t = t_start : param.nTrialsPerBlock
        
        trialNum = trialNum + 1;
        [data exitNow] = sfa_expt2_runTrial(window, param, stim, trialNum, data);
        
        % save progress so far in case we need to restart later
        save([pwd '\data\online_backup\data_sofar.mat'], 'data', 'trialNum', 'blockNum');
        
        if exitNow
            data.timing.blockEnd{blockNum} = clock;
            return; 
        end
    
    end
    data.timing.blockDur(blockNum) = (GetSecs-t0)/60;
    data.timing.blockEnd{blockNum} = clock;
    
    if trialNum < stim.nTrials
        
        if param.isMEG
            blockText = 'Break time!\n\nRelax for a moment.';
        else
            blockText = 'Break time!\n\nRelax for a moment. Press any key to continue.';
        end
        
        if blockNum == 1
            bltext1 = 'block';
        else
            bltext1 = 'blocks';
        end
        
        if stim.nBlocks - blockNum == 1
            bltext2 = 'block';
        else
            bltext2 = 'blocks';
        end
        
        blockText = [blockText '\n\n' num2str(blockNum) ' ' bltext1 ' down, ' num2str(stim.nBlocks - blockNum) ' ' bltext2 ' to go.'];
        DrawFormattedText(window, blockText, sx, sy, param.fontColor);
        Screen('Flip',window);
        
        WaitSecs(2);

        
        if param.isMEG
            % wait indefinitely, to allow experimenter to get MEG computer
            % ready for next recording
            bkey = recordValidKeys(GetSecs, Inf, param.keyboardNumber, {param.readyKey param.exitKey});
            if strcmp(bkey, param.exitKey)
                exitNow = 1;
                return
            end            
            
            blockText = 'Press any key to continue.'; 
            DrawFormattedText(window, blockText, sx, sy, param.fontColor);
            Screen('Flip',window);

            WaitSecs(1);
            KbWait;
        
        else
            bkey = recordValidKeys(GetSecs, Inf, param.keyboardNumber, []);
            if strcmp(bkey, param.exitKey)
                exitNow = 1;
                return
            end            
        
        end
        
        Screen('Flip',window);
        WaitSecs(1);
        
    else
        blockText = 'All done!';
        DrawFormattedText(window, blockText, 'center', 'center', param.fontColor);
        Screen('Flip',window);
        
        WaitSecs(.5);
%         if param.isMEG
%             recordValidKeys(GetSecs, Inf, param.keyboardNumber, param.readyKey);
%         else
%             KbWait;
%         end
        KbWait;
        
    end
end