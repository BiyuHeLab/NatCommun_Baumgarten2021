function NASTD_MEG_SHI_KprimeSummaryStruct...
    (subs, tonedur_text, ...
    paths_NASTD_MEG)
%Aim: Create summary structures containing either experimental or shuffled
%null distribution K-prime values
%Output:
%1)Experimental k-prime vals
%1.1) Single-subject summary file (Kprime_AllSub.Exp({i_tone}{i_win}(i_sub,sensor)) from averaged combined sets
%1.2) Group-Avg summary file (Kprime_GAvg.Exp ({i_tone}{i_win}sensor) = across-subject-averaged experimental k-prime values)

%2) Shuffled k-prime vals
%2.1) Single-subject summary file (kShuffleData_SummarySSubs{i_sub,i_tone}{1,i_win})
%all single-subject shuffled k-prime values, either averaged across sets or both sets appended

%% 0) Specify vars, paths, and setup fieldtrip
addpath(genpath(paths_NASTD_MEG.ScriptsDir));
path_outputdata = [paths_NASTD_MEG.Current_outputdata '/SummaryStruct/'];

%% 1) Create summary files
for i_tonedur = 1:length(tonedur_text)
    %Tone duration condition
    if strcmp(tonedur_text{i_tonedur},'0.15')
        tonedur_title = '0.15sTD';
    elseif  strcmp(tonedur_text{i_tonedur},'0.3')
        tonedur_title = '0.3sTD';
    elseif strcmp(tonedur_text{i_tonedur},'0.6')
        tonedur_title = '0.6sTD';
    end
    
    for i_sub = 1:length(subs)
        %Load in single-subject exp + shuff K-prime vals
        path_inputdata = [paths_NASTD_MEG.Current_outputdata 'ExpK/' tonedur_title '/' subs{i_sub} '/'];
        load([path_inputdata subs{i_sub} '_ExpvsShuffKprime_' tonedur_title '.mat']);
        
        %place single-subject exp + shuff K-prime vals in common struct
        for i_win = 1:length(Kprime.Exp.avgFolds.Kprime.Avg)
            Kprime_AllSub.Exp{i_tonedur}{i_win}(i_sub,:) = Kprime.Exp.avgFolds.Kprime.Avg{i_win};
            Kprime_AllSub.Shuff{i_sub,i_tonedur}{1,i_win} = Kprime.Shuff.avgFolds.NullDist{i_win};
        end
        
        clear Kprime %clean-up
        
    end
    %Compute GroupAvg for EXP and SHUFF data and place it in common struct
    for i_win = 1:length(Kprime_AllSub.Exp{i_tonedur})
        Kprime_GAvg.Exp{i_tonedur}{i_win}(1,:) = mean(Kprime_AllSub.Exp{i_tonedur}{i_win});
    end
end

%% 2) Save Group-level summary output
savefile = [path_outputdata 'KprimeSummaryStruct_SSubsGAvg.mat'];
save(savefile, ...
    'Kprime_GAvg', 'Kprime_AllSub',...
    '-v7.3');
end
