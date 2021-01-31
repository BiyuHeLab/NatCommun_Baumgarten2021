function NASTD_MEG_SHI_CompExpK_perBeta...
    (sub, betalevel_input, tonedur_text, ...
    NumFolds, ...
    paths_NASTD_MEG)

%Aim: Compute k-values with k-fold cross-validation procedure across all
%trials with a specific beta level (across tone durations)

%% 0) specify paths and set up vars
addpath(genpath(paths_NASTD_MEG.ScriptsDir));
NASTD_MEG_SubInfo

%Data input
loadfile_MEG = [si.path_ica  'data_clean.mat'];%path to single-subject MEG data
loadfile_behav = si.path_behav;%path to single-subject behavioral/stimulus data

%Data output
path_save = [paths_NASTD_MEG.Current_outputdata 'ExpK/SeqBeta_' betalevel_input '/'];
mkdir(path_save)

%% 1) load preprocessed MEG data
tic
disp('loading...')
load(loadfile_MEG);
load(loadfile_behav);
disp(['done loading in ' num2str(toc) ' sec'])

inputData_behav = data;
data = [];

for i_tonedur = 1:length(tonedur_text)
    % Loop across tonedur - trial-readout requires TD specification due to
    %TD-specific trial-length. Thus, we first read out, assign to fold, and
    %average over TW for each TD seperately. Then, we append the ERF-results across TD.
    
    data = inputData_behav;
    
    %% 2) Bring data into required format
    %Extract trials from clean MEG data
    data_MEG = extract_trials(data_clean, sub, str2double(tonedur_text{i_tonedur}));
    
    %Filter behavioral data to exclude trials that don't have MEG recording
    numoftrials = num2str(length(data.trialNum)); 
    
    dat_fields  = fieldnames(data);
    for i = 1:length(dat_fields)
        if eval(['length(data.' dat_fields{i} ') == ' numoftrials])
            eval(['data.' dat_fields{i} ' = data.' dat_fields{i} '(si.good_trials);']);
        end
    end
    
    stim_fields = fieldnames(data.stim);
    for i = 1:length(stim_fields)
        if eval(['length(data.stim.' stim_fields{i} ') == ' numoftrials])
            eval(['data.stim.' stim_fields{i} ' = data.stim.' stim_fields{i} '(si.good_trials);']);
        end
    end
    
    %Filter data_MEG for subjects where blocks have to be taken out
    dat_MEG_fields  = fieldnames(data_MEG);
    for i = 1:length(dat_MEG_fields)
        if eval(['length(data_MEG.' dat_MEG_fields{i} ') ==  ' numoftrials])
            eval(['data_MEG.' dat_MEG_fields{i} ' = data_MEG.' dat_MEG_fields{i} '(si.good_trials);']);
        end
    end
    
    %Filter behavioral data & MEG data for trials corresponding to
    %a particular tone length (should be 108 trials per tone dur; tone length ID: 0.15 = 1, 0.3 = 2, 0.6 = 3)
    numoftrials = num2str(length(data.trialNum)); %recount number of trials in case some were rejected above
    f = [];
    f = data.stim.toneDur == str2double(tonedur_text{i_tonedur}); %filter to select trials for certain tone dur
    
    dat_fields  = fieldnames(data);
    for i = 1:length(dat_fields)
        if eval(['length(data.' dat_fields{i} ') == ' numoftrials])
            eval(['data.' dat_fields{i} ' = data.' dat_fields{i} '(f);']);
        end
    end
    
    stim_fields = fieldnames(data.stim);
    for i = 1:length(stim_fields)
        if eval(['length(data.stim.' stim_fields{i} ') == ' numoftrials])
            eval(['data.stim.' stim_fields{i} ' = data.stim.' stim_fields{i} '(f);']);
        end
    end
    
    MEG_fields = fieldnames(data_MEG);
    for i = 1:length(MEG_fields)
        if eval(['length(data_MEG.' MEG_fields{i} ') == ' numoftrials])
            eval(['data_MEG.' MEG_fields{i} ' = data_MEG.' MEG_fields{i} '(f);']);
        end
    end
    
    %a particular sequence beta ID
    numoftrials = num2str(length(data.trialNum)); %recount number of trials in case some were rejected above
    f = [];
    f = data.stim.beta == str2double(betalevel_input); %filter to select trials for certain sequence beta
    
    dat_fields  = fieldnames(data);
    for i = 1:length(dat_fields)
        if eval(['length(data.' dat_fields{i} ') == ' numoftrials])
            eval(['data.' dat_fields{i} ' = data.' dat_fields{i} '(f);']);
        end
    end
    
    stim_fields = fieldnames(data.stim);
    for i = 1:length(stim_fields)
        if eval(['length(data.stim.' stim_fields{i} ') == ' numoftrials])
            eval(['data.stim.' stim_fields{i} ' = data.stim.' stim_fields{i} '(f);']);
        end
    end
    
    MEG_fields = fieldnames(data_MEG);
    for i = 1:length(MEG_fields)
        if eval(['length(data_MEG.' MEG_fields{i} ') == ' numoftrials])
            eval(['data_MEG.' MEG_fields{i} ' = data_MEG.' MEG_fields{i} '(f);']);
        end
    end
    
    %Check MEG data for NaNs and remove those trials
    numoftrials = num2str(length(data.trialNum));
    NaNtrial = [];
    for i_trial = 1:length(data_MEG.trial)
        NaNfind_dataMEG{i_trial} = find(isnan(data_MEG.trial{i_trial}));
        if ~isempty(NaNfind_dataMEG{i_trial})
            NaNtrial(1,:) = i_trial;
        end
    end
    
    if ~isempty(NaNtrial)
        disp(['NaN entries found - removing corrupted trials [' num2str(NaNtrial) ']']);
        f = ~(1:length(data.trialNum) == NaNtrial); %filter to select trials for certain beta level

        dat_fields  = fieldnames(data);
        for i = 1:length(dat_fields)
            if eval(['length(data.' dat_fields{i} ') ==' numoftrials]) 
                eval(['data.' dat_fields{i} ' = data.' dat_fields{i} '(f);']);
            end
        end

        stim_fields = fieldnames(data.stim);
        for i = 1:length(stim_fields)
            if eval(['length(data.stim.' stim_fields{i} ') ==' numoftrials]) 
                eval(['data.stim.' stim_fields{i} ' = data.stim.' stim_fields{i} '(f);']);
            end
        end

        MEG_fields = fieldnames(data_MEG);
        for i = 1:length(MEG_fields)
            if eval(['length(data_MEG.' MEG_fields{i} ') ==' numoftrials])
                eval(['data_MEG.' MEG_fields{i} ' = data_MEG.' MEG_fields{i} '(f);']);
            end
        end
    end        
    
    %% 3) Define folds & allocate trials to folds
    %Determine overall trial number
    NumTrials = length(data_MEG.trial);
    
    %Determine number of folds based on trial number (should be >10,
    %otherwise reduce number of folds)
    if (NumTrials*3)/NumFolds < 10 %*3 here because subsequently, we will combine TD
        NumFolds = NumFolds-1;
        disp('-- Number of folds changed due to insufficient number of reps per uSeq --')
    end
    disp([num2str(NumFolds) ' folds used for cross-validation'])
    
    %Determine how many overall trials are going to be placed in each fold
    AvgTrials_perFold = NumTrials/NumFolds;
    Remaining_trials = NumTrials;
    for i_fold = 1:NumFolds
        if Remaining_trials >= Remaining_trials - (Remaining_trials - ceil(AvgTrials_perFold))
            Trials_perFold(1, i_fold) = Remaining_trials - (Remaining_trials - ceil(AvgTrials_perFold));
            Remaining_trials =  Remaining_trials - ceil(AvgTrials_perFold);
        else
            Trials_perFold(1, i_fold) = Remaining_trials;
        end
    end
    disp([num2str(Trials_perFold) ' trials per fold per TD used for cross-validation'])
    
    %Shuffle trial order within each SeqBeta/TD combination
    trial_index_orig = 1:NumTrials;
    trial_index_shuffled = Shuffle(trial_index_orig,1);
    
    %Allocate shuffled trial to folds
    counter = 1;
    for i_fold = 1:NumFolds
        index_Trials_perFold{i_fold} = trial_index_shuffled(counter:(counter + Trials_perFold(i_fold) -1));
        %Index relative to SeqBeta-TD trials (36)
        counter = counter + Trials_perFold(i_fold);
    end
    
    %Store selected reps in summary file
    for i_fold = 1:NumFolds
        Info_TrialSelection.NumberRepsperFold.PerTD(i_tonedur,:) = Trials_perFold;
        
        Info_TrialSelection.TrialIndexperFold_SeqBetaTDDataSet.PerTD{i_tonedur}{i_fold} = ...
            index_Trials_perFold{i_fold};%Index relative to SeqBeta-TD trials (36 trials)
        Info_TrialSelection.TrialIndexperFold_WholeDataSet.PerTD{i_tonedur}{i_fold} = ...
            data.trialNum(index_Trials_perFold{i_fold});%Index relative to whole data set (324 trials)
        %necessary to get tone sequence data for fold-specific trial later
    end
    index_Trials_perFold = [];
    
    %% 4) For all TD-specific trials, average neural activity within each TW.
    %This results in a struct that can be appended across TD. After
    %that, looping across folds and test/train-trial selection can start.
    
    %Determine window-length parameters and tone start/endings
    nTrials = length(data_MEG.trial);
    nSensors = size( data_MEG.trial{1}, 1 );
    
    sampleFreq      = data_MEG.fsample;
    toneDur_inSecs  = str2num(tonedur_text{i_tonedur});
    nSamplesPerTone = toneDur_inSecs * sampleFreq;
    
    nSamplesPerSeries = length(data_MEG.trial{1}(1,:));
    
    win_size    = 30;
    win_overlap = 0;
    
    t_start_ind(1) = find(data_MEG.time{...
        Info_TrialSelection.TrialIndexperFold_SeqBetaTDDataSet.PerTD...
        {i_tonedur}{1}(1)} == 0); %Any TD-specific trial
    t_end_ind(1)   = t_start_ind(1) + nSamplesPerTone  - 1;
  
    for i_tone = 2:34 %loop across tones
        t_start_ind(i_tone) = t_end_ind(i_tone-1) + 1; %Start index is end index of previous tone +1
        t_end_ind(i_tone)   = t_start_ind(i_tone) + nSamplesPerTone - 1;
        %End index/final data point of i_tone is start+number of samples per tone -1
    end
    
    windows = [1 win_size];
    while windows(end,end) < nSamplesPerTone
        windows = [windows; windows(end,:) + (win_size - win_overlap)];
    end
    
    if windows(end,end) > nSamplesPerTone
        windows(end,:) = [];
    end
    
    %Compute window start/end time in ms for each time window
    windows_inms = (windows / sampleFreq) * 1000;
    numWindows(1,i_tonedur) = length(windows);
    
    min_numWindows = min(numWindows);
    
    %Compute across-time-point-averaged MEG data for each sensor, time window, trial, and tone
    for toneIndex = 1:33 %loop across tones
        for i_trial = 1:nTrials %loop across training trials
            % get MEG data (sensors*time points) for current tone
            dataMEG_perTone_perTrial = data_MEG.trial{i_trial}(:, t_start_ind(toneIndex) : t_end_ind(toneIndex));
            
            % get mean ERF activity for all windows for current tone
            for i_win = 1:min_numWindows
                ind_start = windows(i_win, 1);
                ind_end   = windows(i_win, 2);
                ERF_perTW{i_tonedur}{i_trial}{i_win}{toneIndex}(:,1) = ...
                    nanmean(dataMEG_perTone_perTrial(:, ind_start:ind_end),2);
                %Output: Within-TW-averaged ERF activity per sensor for a specific seqbeta,
                %each TD, each Trial (all folds), each window (1:3),each tone (1:33)
            end
            CorrespondingTrialIndex(i_tonedur,i_trial) = data.trialNum(i_trial);
            CorrespondingToneSeqInformation{i_tonedur}{i_trial} = data.stim.series_f{i_trial};
        end
    end
end

%% 5) Trial-to-fold-assignment
%Combine Trial-to-fold-assignment
Info_TrialSelection.NumberRepsperFold.Total = sum(Info_TrialSelection.NumberRepsperFold.PerTD,1);

%Specify trial-to fold assignment for each fold-loop
Info_TrialSelection.TrialIndexperSeqBeta_WholeDataSet = [];
for i_fold = 1:NumFolds
    
    Info_TrialSelection.TrialIndexperFold_WholeDataSet.AcrossTD{i_fold} = [];
    
    for i_tonedur = 1:length(tonedur_text)
        %Combine trial indices across TD
        Info_TrialSelection.TrialIndexperFold_WholeDataSet.AcrossTD{i_fold} = ...
            [Info_TrialSelection.TrialIndexperFold_WholeDataSet.AcrossTD{i_fold}, ...
            Info_TrialSelection.TrialIndexperFold_WholeDataSet.PerTD{i_tonedur}{i_fold}];
        Info_TrialSelection.TrialIndexperSeqBeta_WholeDataSet = ...
            [Info_TrialSelection.TrialIndexperSeqBeta_WholeDataSet, ...
            Info_TrialSelection.TrialIndexperFold_WholeDataSet.PerTD{i_tonedur}{i_fold}];
    end
end

%Allocate trials to train vs. test sets
for i_fold = 1:NumFolds
    disp(['Current Test Fold = ' num2str(i_fold) '/' num2str(NumFolds)])
    
    %Determine trial indices (whole data set) for train/test sets
    index_test = Info_TrialSelection.TrialIndexperFold_WholeDataSet.AcrossTD{i_fold}; %test data = current fold (i_fold)
    index_train = setdiff(Info_TrialSelection.TrialIndexperSeqBeta_WholeDataSet, ...
        Info_TrialSelection.TrialIndexperFold_WholeDataSet.AcrossTD{i_fold},'stable'); %train data = all other folds
    
    %Select respective trials (MEG and ToneSeqInformation) and place them in common test/train struct
    %Test
    trial_counter = 0;
    for i_trial = index_test
        trial_counter = trial_counter+1;
        [row,col] = find(CorrespondingTrialIndex==i_trial); %read out TD(row) and trialnum(column) of current trial
        
        for i_win = 1:min_numWindows
            for toneIndex = 1:33 %loop across tones
                
                erf_win_test{i_win}{toneIndex}(:, trial_counter) = ERF_perTW{row}{col}{i_win}{toneIndex}(:,1);
            end
        end
        tone_pitch_test_orig(trial_counter, :) = log( inputData_behav.stim.series_f{i_trial} );
    end
    tone_pitch_test = tone_pitch_test_orig(:, [1:end-2]);
    
    %Train
    trial_counter = 0;
    for i_trial = index_train
        trial_counter = trial_counter+1;
        [row,col] = find(CorrespondingTrialIndex==i_trial); %read out TD(row) and trialnum(column) of current trial
        
        for i_win = 1:min_numWindows
            for toneIndex = 1:33 %loop across tones
                
                erf_win_train{i_win}{toneIndex}(:, trial_counter) = ERF_perTW{row}{col}{i_win}{toneIndex}(:,1);
            end
        end
        tone_pitch_train_orig(trial_counter, :) = log( inputData_behav.stim.series_f{i_trial} );
    end
    tone_pitch_train = tone_pitch_train_orig(:, [1:end-2]);
    
    %% 6) Compute k values based on traning set
    %compute k as linear regression (across trials) between sensor-wise, window-wise ERF and
    %tone pitch for tones 32:16 (17 tones total)
    ks = 1:16; %1:16 regression terms/k-values/models - all potential k vals
    
    for i_k = 1:length(ks) %index/current k val (i.e., amount of k models)
        disp(['SeqBeta: ' betalevel_input ' - All TD combined - Training Set (Fold: ' num2str(setdiff(1:NumFolds,i_fold)) ') - Current k-model: ' num2str(i_k)])
        for i_win = 1:min_numWindows
            for i_sensor = 1:nSensors
                
                clear dv iv_cell
                dv = [];
                for toneIndex = 16:32 %17 entries
                    dv = [dv; erf_win_train{i_win}{toneIndex}(i_sensor, :)']; %Column vector with MEG activity for time window at current tone
                    %Output: Column with sensor-specific, time-window averaged MEG activity/ERF for all trials for specific window and toneIndex
                    
                    if toneIndex == 16 %for tone index 16 (i.e, first index entry), start new row
                        for kk = 1:ks(i_k)%loop over k values from 1 to current k val
                            iv_cell{kk} = tone_pitch_train(:, toneIndex - kk + 1); %take tone pitch for all trials for curent tone and also for kk tones back
                        end
                        
                    else %for all other tone indices, append above matrix
                        for kk = 1:ks(i_k)
                            iv_cell{kk} = [iv_cell{kk}; tone_pitch_train(:, toneIndex - kk + 1)];
                        end
                    end
                end
                
                iv = [];
                for kk = 1:ks(i_k)
                    iv(:, kk) = iv_cell{kk}; %conversion from cell to matrix
                    %Output: matrix with tone pitch frequencies;
                    %rows = appended trials per tone (tone 16 to 32),
                    %columns = model order (i.e., tones back (0 to 15 previous
                    %tones)
                end
                                
                % linear regression
                %i.e., for each i_k round/model order, an additional previously presented tone is added
                %and linear regression is computed
                ss = regstats(dv, iv, 'linear',  {'tstat', 'fstat', 'r', 'rsquare', 'adjrsquare'});
                
                r = ss.r; %residuals
                n = length(r);
                sig2_train = sum(r.^2) / n; %sum of squared residuals divided by length of model order
                ss = rmfield(ss, 'r');
                ss.sig2_train = sig2_train;
                
                stats_linear{i_fold}{i_win, i_k, i_sensor} = ss;
                
                clear ss
            end            
        end
    end
    
    %% 7) Characterize model fit on test set
    %compute k as linear regression (across trials) between sensor-wise, window-wise ERF and
    %tone pitch for tones 32:16 (17 tones total), similar to above
    for i_k = 1:length(ks)
        disp(['SeqBeta: ' betalevel_input ' - All TD combined - Test Set (Fold: ' num2str(i_fold) ') - Current k-model: ' num2str(i_k)])
        for i_win = 1:min_numWindows
            for i_sensor = 1:nSensors
                
                clear dv iv_cell
                dv = [];
                for toneIndex = 16:32
                    dv = [dv; erf_win_test{i_win}{toneIndex}(i_sensor, :)'];
                    
                    if toneIndex == 16
                        for kk = 1:ks(i_k)
                            iv_cell{kk} = tone_pitch_test(:, toneIndex - kk +1);
                        end
                        
                    else
                        for kk = 1:ks(i_k)
                            iv_cell{kk} = [iv_cell{kk}; tone_pitch_test(:, toneIndex - kk +1)];
                        end
                    end
                end
                
                iv = [];
                for kk = 1:ks(i_k)
                    iv(:, kk) = iv_cell{kk};
                end
                
                %7.2 Find error of training set on test set
                %Take stats from training set (for current time window, model order, and sensor)
                %and apply them to test set
                ss = stats_linear{i_fold}{i_win, i_k, i_sensor};
                
                beta = ss.tstat.beta';
                
                try
                    for jj = 1:length(dv) %for each trial across all tones from 16:32
                        %compute estimated ERF for each trial across all tones from 16:32
                        %based on the beta weights*tone pitches
                        dv_est(jj) = beta(1) + sum( beta(2:end) .* iv(jj,:) );
                    end
                    
                catch
                    keyboard
                end
                r = dv - dv_est';%compute residuals/difference between estimated and real test set ERF activity for each trial across all tones from 16:32
                
                n = length(r);
                sig2_test = sum(r.^2) / n;
                ss.sig2_test = sig2_test;
                
                % linear regression on test set
                ss2 = regstats(dv, iv, 'linear',  {'tstat', 'fstat', 'r', 'rsquare', 'adjrsquare'});
                
                r = ss2.r;
                n = length(r);
                sig2_test_regr = sum(r.^2) / n;
                
                ss.sig2_test_regr = sig2_test_regr;
                
                stats_linear{i_fold}{i_win, i_k, i_sensor} = ss; %Add test stats to common stats output var
                
                clear ss
            end            
        end
    end
    
end %End Fold_loop

%% 8) Determine Kprime as 'winning' k-model (model with minimal sum squared residuals)
%8.1 Restructure regstats data and read out relevant stats for model selection from even set
n_k      = size(stats_linear{1},2); %number model order;
k_list = 1:n_k;

for i_win = 1:numWindows
    for i_sensor = 1:nSensors
        for i_fold = 1:NumFolds
            
            for i_k = 1:n_k
                sumsquaredRes_proxy(i_k) = stats_linear{i_fold}{i_win, i_k, i_sensor}.sig2_test;
                %sum of squared residuals divided by length of model order for each k (per sensor,time win/fold)
            end
            
            %8.2 Compute criteria for model selection per set and 
            %place timewin/sensor wise parameter from all folds in common struct            
            [min_sumsqaredRes, index_min_sumsquaredRes] = min( sumsquaredRes_proxy );
            %read out minimum squared residuals divided by length of model
            %order and the respective position in array per sensor/time win/fold
            
            CombFolds_kprime_allFolds{i_win}(i_fold,i_sensor) = k_list(index_min_sumsquaredRes) - 1;
            %k-prime-value (i.e., preferred model order/Number of tones
            %back in sequence history best explaining MEG data;
            %corrected by -1 to get kprime value (0-15), not index/position (1-16))
            
            CombFolds_minsumsquaredRes_allFolds{i_win}(i_fold,i_sensor) = min_sumsqaredRes;
            %minimum squared residuals divided by length of model order of K-value (i.e., preferred model order/Number of tones back)
            
            CombFolds_sumsquaredRes_allK_allFolds{i_win}{i_sensor}(i_fold, :) = ...
                sumsquaredRes_proxy;
            %sum of squared residuals divided by length of model order,
            %stored in common matrix across folds
            
            %8.3 Average beta regression weights for each k-model across sets
            %note: there are k+1 beta weights for each k-model, with
            %beta weight (1) being the offset term. We interested in
            %beta weightes (2+) since they resemble the beta weights
            %per k-model
            for i_k = 1:n_k
                CombFolds_BetaRegWeights_allFolds{i_win, i_k, i_sensor}(i_fold,:) =  ...
                    stats_linear{i_fold}{i_win, i_k, i_sensor}.tstat.beta';
            end
            
            %12.4 Clean-up
            clear sumsquaredRes_proxy min_sumsqaredRes index_min_sumsquaredRes
            
        end
    end
end


%% 9) Combine Kprime across folds by averaging values across cross-validation folds
%i.e., average model selection criteria by averaging across folds
for i_win = 1:numWindows
    for i_sensor = 1:nSensors
        
        %9.1 Average and compute STD across folds (for each sensor/time win)
        %Kprime
        CombFolds_kprime_AVGacrossFolds{i_win}(i_sensor) = mean(CombFolds_kprime_allFolds{i_win}(:,i_sensor));
        CombFolds_kprime_STDacrossFolds{i_win}(i_sensor) = std(CombFolds_kprime_allFolds{i_win}(:,i_sensor));
        %SSR Kprime
        CombFolds_minsumsquaredRes_AVGacrossFolds{i_win}(i_sensor) = mean(CombFolds_minsumsquaredRes_allFolds{i_win}(:,i_sensor));
        CombFolds_minsumsquaredRes_STDacrossFolds{i_win}(i_sensor) = std(CombFolds_minsumsquaredRes_allFolds{i_win}(:,i_sensor));
        %SSR all K
        CombFolds_sumsquaredRes_allK_AVGacrossFolds{i_win}(i_sensor, :) = mean(CombFolds_sumsquaredRes_allK_allFolds{i_win}{i_sensor});
        CombFolds_sumsquaredRes_allK_STDacrossFolds{i_win}(i_sensor, :) = std(CombFolds_sumsquaredRes_allK_allFolds{i_win}{i_sensor});
        
        %Beta Regression Weights
        %all K
        for i_k = 1:n_k
            CombFolds_BetaRegWeights_AVGacrossFolds{i_win, i_k, i_sensor} = mean(CombFolds_BetaRegWeights_allFolds{i_win, i_k, i_sensor});
            CombFolds_BetaRegWeights_STDacrossFolds{i_win, i_k, i_sensor} = std(CombFolds_BetaRegWeights_allFolds{i_win, i_k, i_sensor});
        end
        
        %K-prime
        CombFolds_kprime_BetaRegWeight_AVGacrossFolds{i_win}{i_sensor} = mean(CombFolds_BetaRegWeights_AVGacrossFolds{i_win, round(CombFolds_kprime_AVGacrossFolds{i_win}(i_sensor) + 1), i_sensor});
        CombFolds_kprime_BetaRegWeight_STDacrossFolds{i_win}{i_sensor} = std(CombFolds_BetaRegWeights_AVGacrossFolds{i_win, round(CombFolds_kprime_AVGacrossFolds{i_win}(i_sensor) + 1), i_sensor});
        %Note: Beta Regression Weights are selected based on 1)
        %across-fold-averaged beta weights and 2) across-fold-averaged Kprime
        %values (here -1 lin reg offset term correction is undone, since we
        %need position in 1:16 aray, not Kprime value)
    end
end

%9.2 Create summary file containing all relevant vars
Kprime_data.perFold.Kprime = CombFolds_kprime_allFolds;
Kprime_data.perFold.SSR_perK = CombFolds_sumsquaredRes_allK_allFolds;
Kprime_data.perFold.SSR_Kprime = CombFolds_minsumsquaredRes_allFolds;
Kprime_data.perFold.BetaRegWeights_perK = CombFolds_BetaRegWeights_allFolds;

Kprime_data.avgFolds.Kprime.Avg = CombFolds_kprime_AVGacrossFolds;
Kprime_data.avgFolds.Kprime.STD = CombFolds_kprime_STDacrossFolds;
Kprime_data.avgFolds.SSR_Kprime.Avg = CombFolds_minsumsquaredRes_AVGacrossFolds;
Kprime_data.avgFolds.SSR_Kprime.STD = CombFolds_minsumsquaredRes_STDacrossFolds;
Kprime_data.avgFolds.SSR_perK.Avg = CombFolds_sumsquaredRes_allK_AVGacrossFolds;
Kprime_data.avgFolds.SSR_perK.STD = CombFolds_sumsquaredRes_allK_STDacrossFolds;
Kprime_data.avgFolds.BetaRegWeights_perK.Avg = CombFolds_BetaRegWeights_AVGacrossFolds;
Kprime_data.avgFolds.BetaRegWeights_perK.STD = CombFolds_BetaRegWeights_STDacrossFolds;
Kprime_data.avgFolds.BetaRegWeights_Kprime.Avg = CombFolds_kprime_BetaRegWeight_AVGacrossFolds;
Kprime_data.avgFolds.BetaRegWeights_Kprime.STD = CombFolds_kprime_BetaRegWeight_STDacrossFolds;

clear CombFolds_*

%% 10) Save output variables
savefile = [path_save sub '_ExpKprime_SeqBeta' betalevel_input '.mat'];
save(savefile, 'Info_TrialSelection','Kprime_data'); %without fold-wise output

end