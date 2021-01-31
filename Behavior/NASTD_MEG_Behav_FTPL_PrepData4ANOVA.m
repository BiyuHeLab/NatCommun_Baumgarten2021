function [subPred, subDurPred] = NASTD_MEG_Behav_FTPL_PrepData4ANOVA(subs, tonedur_text, paths_NASTD_MEG)

%Reads out FTPLrating data looped across subjects and restructures it in a
%way that anova functions can use (i.e., FTPL ratings ordered by p*34, 34,
%tone dur), both across tone dur (subPred) and per tone dur (subDurPred)

possible_predp34 = {'low' 'med' 'high'}; %math. expected final tone pitches
possible_p34 = {'-3' '-2' '-1' '1' '2' '3'}; %actually presented final tone pitches
possible_beta = {'0.5' '0.99' '1.5'}; %sequence beta

subPred = []; %empty placeholder for subject responses

for i_sub = 1:length(subs)
    
    sub = subs{i_sub}
    
    %Option1: Behavioral data from training session for sub 1-7
    load([paths_NASTD_MEG.Rawdata.Behav 'subjects/sfa_expt3_' sub '.mat']);
 
    %     %Option2: Behavioral data from MEG recording day
    %     NASTD_MEG_subjectinfo %load subject info file (var: si)
    %     loadfile_behav = si.path_behav; %path to single-subject behavioral/stimulus data (MEG recording date)
    %     load(loadfile_behav);
    
    vecPredRating = []; %FTPL rating
    vec_tonedur = []; %Tone duration
    vec_predp34 = []; %math. expected final tone (-1 = low, 0 = medium, 1 = high)
    vec_p34 = []; %actually presented final tone frequency
    vec_beta = []; %sequence beta value
    counter_run = 0;
    
    for i_tonedur = 2:length(tonedur_text) %no all condition
        toneDur = tonedur_text{i_tonedur};
        
        for i_expectedFTP = 1:length(possible_predp34)
            expected = possible_predp34{i_expectedFTP};
            
            for i_actualFTP = 1:length(possible_p34)
                actual = possible_p34{i_actualFTP};
                
                for i_beta = 1:length(possible_beta)
                    beta = possible_beta{i_beta};
                    
                    counter_run = counter_run + 1;                    
                    
                    switch toneDur %filter tone duration
                        case '0.15'
                            filter_toneDur = stim.toneDur == 0.15;                            
                        case '0.3'
                            filter_toneDur = stim.toneDur == 0.3;                            
                        case '0.6'
                            filter_toneDur = stim.toneDur == 0.6;                            
                        otherwise
                            filter_toneDur = ones(1,length(data.trialNum));
                    end
                    
                    
                    switch expected %filter math. expected final tone
                        case 'low'
                            filter_expected = stim.predID == -1;
                        case 'med'
                            filter_expected = stim.predID == 0;
                        case 'high'
                            filter_expected = stim.predID == 1;
                    end
                    
                    switch actual %filter actually presented final tone
                        case '-3'
                            filter_actual = stim.finalID == -3;
                        case '-2'
                            filter_actual = stim.finalID == -2;
                        case '-1'
                            filter_actual = stim.finalID == -1;
                        case '1'
                            filter_actual = stim.finalID == 1;
                        case '2'
                            filter_actual = stim.finalID == 2;
                        case '3'
                            filter_actual = stim.finalID == 3;
                    end
                    
                    switch beta %filter beta
                        case '0.5'
                            filter_beta = stim.betaID == 1;
                        case '0.99'
                            filter_beta = stim.betaID == 2;
                        case '1.5'
                            filter_beta = stim.betaID == 3;
                    end
                    
                    filter = filter_toneDur & filter_expected & filter_actual & filter_beta;
                    
                    %% 2) Place trial selection in new structure used for MATLAB ANOVA
                    vecPredRating = [vecPredRating; data.resp_prob(filter ==1)'];
                    %vector reading out the subject's final tone likelihood rating for each trial in
                    %order determined by above loops (i.e., hierarchical order)
                    
                    nTrialsPerCond = length(data.resp_prob(filter ==1));
                    
                    %Update row
                    if counter_run == 1 %if first entry, start on top
                        rowCount = 1;
                    else
                        rowCount = lastRow + 1; %otherwise append
                    end
                    lastRow = rowCount + nTrialsPerCond - 1;
                    
                    %Copy order of loop vars into vectors (i.e., hierarchical sorting for
                    %1) tone duration, 2) math. expected final tone, 3) actually presented
                    %final tone
                    [vec_tonedur{rowCount: lastRow}] = deal(toneDur);
                    [vec_predp34{rowCount: lastRow}] = deal(expected);
                    [vec_p34{rowCount: lastRow}] = deal(actual);
                    [vec_beta{rowCount: lastRow}] = deal(beta);
                    
                end
            end
        end
    end
    
    %Output: single-subject cell containing vectors with all trials,
    %hiearchically sorted by 1) tone dur (1:3), 2) Expected tone (1:3), 3)
    %Actual tone (1:6), 4) sequence beta
    subPred{i_sub}.PredRating = vecPredRating;
    subPred{i_sub}.tonedur = vec_tonedur';
    subPred{i_sub}.predp34 = vec_predp34';
    subPred{i_sub}.p34 = vec_p34';
    subPred{i_sub}.beta = vec_beta';
    
    %Replace -1 entires (missed responses) with NaNs
    i_missedResp = find(subPred{i_sub}.PredRating == -1);
    subPred{i_sub}.PredRating(i_missedResp) = NaN;
    
    %Subsitute numerical values
    for i_trial = 1:length(subPred{i_sub}.tonedur)
        if strcmp(subPred{i_sub}.tonedur(i_trial),'0.15')
            subPred{i_sub}.tonedur_numerical(i_trial,1) = 1;
        elseif  strcmp(subPred{i_sub}.tonedur(i_trial),'0.3')
            subPred{i_sub}.tonedur_numerical(i_trial,1) = 2;
        elseif  strcmp(subPred{i_sub}.tonedur(i_trial),'0.6')
            subPred{i_sub}.tonedur_numerical(i_trial,1) = 3;
        end
        
        if strcmp(subPred{i_sub}.predp34(i_trial),'low')
            subPred{i_sub}.predp34_numerical(i_trial,1) = 1;
        elseif  strcmp(subPred{i_sub}.predp34(i_trial),'med')
            subPred{i_sub}.predp34_numerical(i_trial,1) = 2;
        elseif  strcmp(subPred{i_sub}.predp34(i_trial),'high')
            subPred{i_sub}.predp34_numerical(i_trial,1) = 3;
        end
        
        if strcmp(subPred{i_sub}.p34(i_trial),'-3')
            subPred{i_sub}.p34_numerical(i_trial,1) = 1;
        elseif  strcmp(subPred{i_sub}.p34(i_trial),'-2')
            subPred{i_sub}.p34_numerical(i_trial,1) = 2;
        elseif  strcmp(subPred{i_sub}.p34(i_trial),'-1')
            subPred{i_sub}.p34_numerical(i_trial,1) = 3;
        elseif strcmp(subPred{i_sub}.p34(i_trial),'1')
            subPred{i_sub}.p34_numerical(i_trial,1) = 4;
        elseif  strcmp(subPred{i_sub}.p34(i_trial),'2')
            subPred{i_sub}.p34_numerical(i_trial,1) = 5;
        elseif  strcmp(subPred{i_sub}.p34(i_trial),'3')
            subPred{i_sub}.p34_numerical(i_trial,1) = 6;
        end
        
        if strcmp(subPred{i_sub}.beta(i_trial),'0.5')
            subPred{i_sub}.beta_numerical(i_trial,1) = 1;
        elseif  strcmp(subPred{i_sub}.beta(i_trial),'0.99')
            subPred{i_sub}.beta_numerical(i_trial,1) = 2;
        elseif  strcmp(subPred{i_sub}.beta(i_trial),'1.5')
            subPred{i_sub}.beta_numerical(i_trial,1) = 3;
        end  
        
    end
end

tone_durs = [0.15 0.3 0.6];

%Seperate data into seperate cells per tone duration
for i_sub = 1:length(subs)
    for i_tone = 1:length(tone_durs)
        
        f_tone = str2num(char(subPred{i_sub}.tonedur)) == tone_durs(i_tone);
        
        subDurPred{i_sub, i_tone}.PredRating    = subPred{i_sub}.PredRating(f_tone);
        subDurPred{i_sub, i_tone}.predp34      = subPred{i_sub}.predp34(f_tone);
        subDurPred{i_sub, i_tone}.p34        = subPred{i_sub}.p34(f_tone);
        subDurPred{i_sub, i_tone}.beta        = subPred{i_sub}.beta(f_tone);
       
    end
end

%Subsitute numerical values
for i_sub = 1:length(subs)
    for i_tone = 1:length(tone_durs)
        for i_trial = 1:length(subDurPred{i_sub, i_tone}.predp34)
            
            if strcmp(subDurPred{i_sub, i_tone}.predp34(i_trial),'low')
                subDurPred{i_sub, i_tone}.predp34_numerical(i_trial,1) = 1;
            elseif  strcmp(subDurPred{i_sub, i_tone}.predp34(i_trial),'med')
                subDurPred{i_sub, i_tone}.predp34_numerical(i_trial,1) = 2;
            elseif  strcmp(subDurPred{i_sub, i_tone}.predp34(i_trial),'high')
                subDurPred{i_sub, i_tone}.predp34_numerical(i_trial,1) = 3;
            end
            
            if strcmp(subDurPred{i_sub, i_tone}.p34(i_trial),'-3')
                subDurPred{i_sub, i_tone}.p34_numerical(i_trial,1) = 1;
            elseif  strcmp(subDurPred{i_sub, i_tone}.p34(i_trial),'-2')
                subDurPred{i_sub, i_tone}.p34_numerical(i_trial,1) = 2;
            elseif  strcmp(subDurPred{i_sub, i_tone}.p34(i_trial),'-1')
                subDurPred{i_sub, i_tone}.p34_numerical(i_trial,1) = 3;
            elseif strcmp(subDurPred{i_sub, i_tone}.p34(i_trial),'1')
                subDurPred{i_sub, i_tone}.p34_numerical(i_trial,1) = 4;
            elseif  strcmp(subDurPred{i_sub, i_tone}.p34(i_trial),'2')
                subDurPred{i_sub, i_tone}.p34_numerical(i_trial,1) = 5;
            elseif  strcmp(subDurPred{i_sub, i_tone}.p34(i_trial),'3')
                subDurPred{i_sub, i_tone}.p34_numerical(i_trial,1) = 6;
            end

            if strcmp(subDurPred{i_sub, i_tone}.beta(i_trial),'0.5')
                subDurPred{i_sub, i_tone}.beta_numerical(i_trial,1) = 1;
            elseif  strcmp(subDurPred{i_sub, i_tone}.beta(i_trial),'0.99')
                subDurPred{i_sub, i_tone}.beta_numerical(i_trial,1) = 2;
            elseif  strcmp(subDurPred{i_sub, i_tone}.beta(i_trial),'1.5')
                subDurPred{i_sub, i_tone}.beta_numerical(i_trial,1) = 3;
            end            
        end
    end
end

end