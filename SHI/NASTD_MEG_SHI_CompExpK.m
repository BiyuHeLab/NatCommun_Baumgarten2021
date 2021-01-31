function NASTD_MEG_SHI_CompExpK...
    (sub, tonedur_text, ...
    NumFolds, ...
    paths_NASTD_MEG)

%Aim: Compute k-values with 6-fold cross-validation procedure as linear
%regression (across trials) between sensor-wise, window-wise MEG
%activity and tone pitch for tones 32:16

%% 0) Specify vars, paths, and setup fieldtrip
addpath(genpath(paths_NASTD_MEG.ScriptsDir));
NASTD_MEG_SubInfo

%Tone duration condition for data load-in
if strcmp(tonedur_text,'0.15')
    tonedur_title = '0.15sTD';
elseif  strcmp(tonedur_text,'0.3')
    tonedur_title = '0.3sTD';
elseif strcmp(tonedur_text,'0.6')
    tonedur_title = '0.6sTD';
end

%Data input
loadfile_MEG = [si.path_ica  'data_clean.mat'];%path to single-subject MEG data
loadfile_behav = si.path_behav;%path to single-subject behavioral/stimulus data

%Data output
path_save = [paths_NASTD_MEG.Current_outputdata 'ExpK/' tonedur_title '/' sub '/'];
mkdir(path_save)

%% 1) load preprocessed MEG data
tic
disp('loading...')
load(loadfile_MEG);
load(loadfile_behav);
disp(['done loading in ' num2str(toc) ' sec'])

%% 2) Bring data into required format
%Extract trials from clean MEG data
data_MEG = extract_trials(data_clean, sub, str2double(tonedur_text));

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

%Determine and label individual sequence ID
%9 unique sequences (3 beta * 3 p*34)* 3 tone dur = 27 trials/blocks, 12
%blocks total = 324 trials total
data.stim.seqID = (1:length(data.stim.series_f));
data.stim.seqID(1:end) = NaN;
seqID_perUniqueSum = [];%Create seqID proxy field


%Read out all trials of first block and determine individual sequence ID
for i_trial = 1:length(data.stim.series_f)/12 %trials per block (27)
    All_trials_block1(i_trial,:) = data.stim.series_f{i_trial}(1:33); %only first 33 tone, since 34th varies
end
seqID_perUniqueSum = unique(sum(All_trials_block1')); %use sum across tone frequency values as metric determinig SeqID
%should result in 9 unique sequences (independent of final tone pitch/p34)

%Individually label all trials based on their across-tone-freq-sum
for i_trial = 1:length(data.stim.series_f)
    for i_seqID = 1:length(seqID_perUniqueSum)
        if sum(data.stim.series_f{i_trial}(1:33)) == seqID_perUniqueSum(i_seqID)
            data.stim.seqID(i_trial) = i_seqID;
        end
    end
end

%Filter trials with selected ToneDur
numoftrials = num2str(length(data.trialNum));
filter_tonedur = data.stim.toneDur == str2double(tonedur_text);

dat_fields  = fieldnames(data);
for i = 1:length(dat_fields)
    if eval(['length(data.' dat_fields{i} ') ==' numoftrials])
        eval(['data.' dat_fields{i} ' = data.' dat_fields{i} '(filter_tonedur);']);
    end
end

stim_fields = fieldnames(data.stim);
for i = 1:length(stim_fields)
    if eval(['length(data.stim.' stim_fields{i} ') ==' numoftrials])
        eval(['data.stim.' stim_fields{i} ' = data.stim.' stim_fields{i} '(filter_tonedur);']);
    end
end

MEG_fields = fieldnames(data_MEG);
for i = 1:length(MEG_fields)
    if eval(['length(data_MEG.' MEG_fields{i} ') ==' numoftrials])
        eval(['data_MEG.' MEG_fields{i} ' = data_MEG.' MEG_fields{i} '(filter_tonedur);']);
    end
end

%Tone Duration Filtered Sanity check: count number of seqID trials 
%(should be 12 per seqID per tone duration condition)
for i_seqID = 1:length(seqID_perUniqueSum)
    if length(find(data.stim.seqID==i_seqID)) ~= 12
        warning0 = ['CAVE: Not all  repeptitions present for uSeq:' ...
            num2str(i_seqID) ' in tonedur condition ' tonedur_title];
        disp(warning0);
    end
end

%Check MEG data for NaNs
NaNtrial = [];
for i_trial = 1:length(data_MEG.trial)
    NaNfind_dataMEG{i_trial} = find(isnan(data_MEG.trial{i_trial}));
    if ~isempty(NaNfind_dataMEG{i_trial})
        NaNtrial(1,:) = i_trial;
    end
end
if ~isempty(NaNtrial)
    newTrialSelection = 1:length(data_MEG.trial);
    newTrialSelection(NaNtrial) = [];
    
    sampleinfo = data_MEG.sampleinfo;
    sampleinfo(NaNtrial) = [];
    %If NaNs present, delete trials with NaNs
    cfg = [];
    cfg.trials =  newTrialSelection;
    data_MEG = ft_selectdata(cfg, data_MEG);
    data_MEG.sampleinfo = sampleinfo;
end

%% 3) Define folds, allocate trials to folds
%Determine uSeq-wise trial selection
Number_uSeq = length(unique(data.stim.seqID));

%Create empty proxies
Number_reps_peruSeq = NaN(1,Number_uSeq);
index_reps_peruSeq = NaN(Number_uSeq,12);
minNumberacrosstonedur_reps_peruSeq = NaN(1,Number_uSeq);

index_reps_peruSeq_noNaN = [];
index_train_peruSeq = [];
index_test_peruSeq = [];
temp1 =  [];

index_train = [];
index_test = [];

%Find trial indices corresponding to reps in each uSeq
for i_seqID = 1:Number_uSeq
    index_reps_peruSeq(i_seqID,1:length(find(data.stim.seqID == i_seqID))) = find(data.stim.seqID == i_seqID);
    
    %Check if any trial (determined by data.stim) is not present in MEG
    for i_entry = 1:length(index_reps_peruSeq(i_seqID,1:length(find(data.stim.seqID == i_seqID))))
        if index_reps_peruSeq(i_seqID,i_entry) > length(data_MEG.trial)
            index_reps_peruSeq(i_seqID,i_entry) = NaN;
        end
    end
    
    %get rid of NaN entries
    counter = 1;
    for i_rep = 1:length(index_reps_peruSeq)
        if ~isnan(index_reps_peruSeq(i_seqID,i_rep))
            index_reps_peruSeq_noNaN{i_seqID}(counter)  = index_reps_peruSeq(i_seqID,i_rep);
            counter = counter+1;
        end
    end
end

%Allocate reps to folds (optimal if identical number of reps per group, but not necessary)
%Update Number of folds based on reps per uSeq (to guarantee that no fold is empty)
for i_seqID = 1:Number_uSeq
    Number_reps_peruSeq(i_seqID) = length(index_reps_peruSeq_noNaN{i_seqID});
end
AvgNumber_reps_peruSeq = ceil(mean(Number_reps_peruSeq));

if AvgNumber_reps_peruSeq < 11 && NumFolds > 2
    NumFolds = NumFolds-1;
    disp('-- Number of folds changed due to insufficient number of reps per uSeq --')
end

disp([num2str(NumFolds) ' folds used for cross-validation'])

%Determine how many reps peruSeq will be placed in each fold
for i_seqID = 1:Number_uSeq
    Reps_perFold{i_seqID} = length(index_reps_peruSeq_noNaN{i_seqID})/NumFolds;
    if mod(Reps_perFold{i_seqID},1) ~= 0
        warning1 = ['Uneven number of reps/fold for uSeq:' num2str(i_seqID)];
        disp(warning1)
    end
end

%Shuffle order of reps before allocating them to folds
for i_seqID = 1:Number_uSeq
    index_reps_peruSeq_noNaN_shuffled{i_seqID} = Shuffle(index_reps_peruSeq_noNaN{i_seqID});
    Info_TrialSelection.Number_Usedreps_peruSeq(1,i_seqID) = length(index_reps_peruSeq_noNaN{i_seqID});
    Info_TrialSelection.TrialIndexUsedreps_ToneDurDataSet_peruSeq{i_seqID} = index_reps_peruSeq_noNaN{i_seqID};
    Info_TrialSelection.TrialIndexUsedreps_WholeDataSet_peruSeq{i_seqID} = data.trialNum(index_reps_peruSeq_noNaN{i_seqID});    
end

%Allocate shuffled reps to folds
%For each fold, a balanced number of uSeqs reps is added. If there are less
%than 12 reps, less are added to the last fold.
counter = 1;
f_unpresenttrials = [];
for i_fold = 1:NumFolds
    index_reps_peruSeq_perfold{i_fold} = [];
    
    for i_seqID = 1:Number_uSeq
        index_selectedtrials = counter:counter+ceil(Reps_perFold{i_seqID}-1);
        for i_selectedtrial = 1:length(index_selectedtrials)
            if index_selectedtrials(i_selectedtrial) > length(index_reps_peruSeq_noNaN_shuffled{i_seqID})
                f_unpresenttrials(i_selectedtrial) = 1;
            end
        end
        f_unpresenttrials = logical(f_unpresenttrials);
        index_selectedtrials(f_unpresenttrials) = [];
        f_unpresenttrials = [];
        
        index_reps_peruSeq_perfold{i_fold} = [index_reps_peruSeq_perfold{i_fold} ...
            index_reps_peruSeq_noNaN_shuffled{i_seqID}(index_selectedtrials)];
    end
    counter = counter + ceil(Reps_perFold{i_seqID});
    
    %Store selected reps in summary file
    Info_TrialSelection.Number_Usedreps_perFold{i_fold} = length(index_reps_peruSeq_perfold{i_fold});
    Info_TrialSelection.TrialIndexUsedreps_ToneDurDataSet_perFold{i_fold} = index_reps_peruSeq_perfold{i_fold};
    Info_TrialSelection.TrialIndexUsedreps_WholeDataSet_perFold{i_fold} = data.trialNum(index_reps_peruSeq_perfold{i_fold});
    
end

%% 4) Loop subsequent computation over folds
for i_fold = 1:NumFolds
    disp(['Current Test Fold = ' num2str(i_fold) '/' num2str(NumFolds)])
    
    %For each fold_loop, specify train and test data
    index_test = index_reps_peruSeq_perfold{i_fold}; %test data = current fold
    index_train = []; %train data = all other folds
    for i_fold2 = 1:NumFolds
        index_train = [index_train index_reps_peruSeq_perfold{i_fold2}];
    end
    for i_rep = 1:length(index_reps_peruSeq_perfold{i_fold})
        index_train(find(index_train == index_reps_peruSeq_perfold{i_fold}(i_rep))) = [];
    end
    
    %Construct corresponding filters
    f_train = [1:length(data_MEG.trial)]*0; %create empty filter-proxies
    f_test = [1:length(data_MEG.trial)]*0;
    
    f_train(index_train) = 1; %fill set-specific chosen trial indices
    f_test(index_test) = 1;
    
    f_train = logical(f_train); %convert arrays to logical
    f_test = logical(f_test);
    
    %% 5) Define training set
    %Copy field struct and input from selected train trials
    data_MEG_train = data_MEG;
    data_MEG_train.time  = data_MEG_train.time(f_train);
    data_MEG_train.trial = data_MEG_train.trial(f_train);
    data_MEG_train.sampleinfo = data_MEG_train.sampleinfo(f_train);
    data_MEG_train.trialinfo  = data_MEG_train.trialinfo(f_train);
    
    data_train = data;
    numoftrials = num2str(length(data.trialNum));
    
    dat_fields  = fieldnames(data_train);
    for i = 1:length(dat_fields)
        if eval(['length(data_train.' dat_fields{i} ') == ' numoftrials])
            eval(['data_train.' dat_fields{i} ' = data_train.' dat_fields{i} '(f_train);']);
        end
    end
    
    stim_fields = fieldnames(data_train.stim);
    for i = 1:length(stim_fields)
        if eval(['length(data_train.stim.' stim_fields{i} ') == ' numoftrials])
            eval(['data_train.stim.' stim_fields{i} ' = data_train.stim.' stim_fields{i} '(f_train);']);
        end
    end
    
    %Copy response/behav data and stim data into MEG struct
    data_MEG_train.behav = data_train;
    data_MEG_train.stim  = data_train.stim;
    
    %% 6) Define test set
    %Copy field struct and input from selected train trials
    data_MEG_test = data_MEG;
    data_MEG_test.time  = data_MEG_test.time(f_test);
    data_MEG_test.trial = data_MEG_test.trial(f_test);
    data_MEG_test.sampleinfo = data_MEG_test.sampleinfo(f_test);
    data_MEG_test.trialinfo  = data_MEG_test.trialinfo(f_test);
    
    data_test = data;
    numoftrials = num2str(length(data.trialNum));
    
    dat_fields  = fieldnames(data_test);
    for i = 1:length(dat_fields)
        if eval(['length(data_test.' dat_fields{i} ') == ' numoftrials])
            eval(['data_test.' dat_fields{i} ' = data_test.' dat_fields{i} '(f_test);']);
        end
    end
    
    stim_fields = fieldnames(data_test.stim);
    for i = 1:length(stim_fields)
        if eval(['length(data_test.stim.' stim_fields{i} ') == ' numoftrials])
            eval(['data_test.stim.' stim_fields{i} ' = data_test.stim.' stim_fields{i} '(f_test);']);
        end
    end
    
    %Copy response/behav data and stim data into MEG struct
    data_MEG_test.behav = data_test;
    data_MEG_test.stim  = data_test.stim;
    
    %% 7) Analyse time window defined MEG data
    %Define analysis parameters
    nTrials_train = length(data_MEG_train.trial);
    nTrials_test  = length(data_MEG_test.trial);
    nSensors      = size( data_MEG_train.trial{1}, 1 );
    
    samplingFreq    = data_MEG_train.fsample;
    toneDur_inSecs  = str2num(tonedur_text);
    nSamplesPerTone = toneDur_inSecs * samplingFreq;
    
    nSamplesPerSeries = length(data_MEG_train.trial{1}(1,:));
    
    %Define time windows used for analysis
    win_size    = 30; 
    win_overlap = 0;
    
    %Define starting point (i.e., first tone)
    t_start_ind(1) = find(data_MEG_train.time{1} == 0);
    t_end_ind(1)   = t_start_ind(1) + nSamplesPerTone -1;
    
    %Read out start and stop samples for each tone (except the first)
    for i_tone = 2:34 %loop across tones
        t_start_ind(i_tone) = t_end_ind(i_tone-1) + 1; %Start index is end index of previous tone +1
        t_end_ind(i_tone)   = t_start_ind(i_tone) + nSamplesPerTone - 1;
        %End index/final data point of i_tone is start+number of samples per tone -1
    end
    
    %Define number, start and end sample of window per tone
    windows = [1 win_size];
    while windows(end,end) < nSamplesPerTone
        windows = [windows; windows(end,:) + (win_size - win_overlap)];
    end
    
    if windows(end,end) > nSamplesPerTone
        windows(end,:) = [];
    end
    
    %Compute window start/end time in ms for each time window
    windows_inms = (windows / samplingFreq) * 1000;
    
    %% 8) Compute across-time-point-averaged MEG data for each sensor, time window, trial, and tone
    %Training set
    for toneIndex = 1:33
        for i_trial = 1:nTrials_train
            % get MEG data (sensors*time points) for current tone
            meg_tone_trial = data_MEG_train.trial{i_trial}(:, t_start_ind(toneIndex) : t_end_ind(toneIndex));
            
            % get mean ERF activity for all windows for current tone
            for i_win = 1:size(windows,1)
                ind_start = windows(i_win, 1);
                ind_end   = windows(i_win, 2);
                
                erf_win_train{i_win}{toneIndex}(:, i_trial) = ...
                    mean(meg_tone_trial(:, ind_start:ind_end),2);
            end
        end
    end
    
    %Test set
    for toneIndex = 1:33
        for i_trial = 1:nTrials_test
            % get MEG data for current tone at current trial
            meg_tone_trial = data_MEG_test.trial{i_trial}(:, t_start_ind(toneIndex) : t_end_ind(toneIndex));
            
            % get mean ERF activity for all windows
            for i_win = 1:size(windows,1)
                ind_start = windows(i_win, 1);
                ind_end   = windows(i_win, 2);
                
                erf_win_test{i_win}{toneIndex}(:, i_trial) = ...
                    mean( meg_tone_trial(:, ind_start:ind_end), 2 );
            end
        end
    end
    
    %% 9) Read out previously-presented tone pitches (i.e., sequence history)
    %read out presented tone pitches for training and test set from stim data
    for i_trial = 1:nTrials_train
        tone_pitch_train_o(i_trial, :) = log( data_MEG_train.stim.series_f{i_trial} );
    end
    
    for i_trial = 1:nTrials_test
        tone_pitch_test_o(i_trial, :)  = log( data_MEG_test.stim.series_f{i_trial} );
    end
    
    tone_pitch_train = tone_pitch_train_o(:, [1:end-2]);
    %take original matrix, delete last two column/tone, because we only use first 32 tone pitches
    tone_pitch_test = tone_pitch_test_o(:, [1:end-2]);
    
    %% 10) Compute k values based on traning set
    %compute k as linear regression (across trials) between sensor-wise,
    %window-wise MEG activity and tone pitch for tones 32:16
    ks = 1:16; %regression terms/k-values/models - all potential k vals
    
    for i_k = 1:length(ks)
        disp(['ToneDur: ' tonedur_title '0ms - Training Set (Fold: ' ...
            num2str(i_fold) ') - Current k-model: ' num2str(i_k)])
        for i_win = 1:size(windows,1)
            for i_sensor = 1:nSensors
                
                clear dv iv_cell
                dv = [];
                for toneIndex = 16:32
                    dv = [dv; erf_win_train{i_win}{toneIndex}(i_sensor, :)'];
                    %Column vector with MEG activity for time window at current tone across training trials
                    %Output: Column with sensor-specific, time-window averaged MEG activity for all trials for specific window and toneIndex
                    
                    if toneIndex == 16
                        for kk = 1:ks(i_k)
                            iv_cell{kk} = tone_pitch_train(:, toneIndex - kk + 1);
                            %take tone pitch for all trials for curent tone and also for kk tones back
                        end
                    else
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
                
                %Compute linear regression
                %i.e., for each i_k round/model order, an additional
                %previously presented tone is added and linear regression is computed
                ss = regstats(dv, iv, 'linear',  {'tstat', 'fstat', 'r', 'rsquare', 'adjrsquare'});
                r = ss.r; %residuals
                n = length(r);
                sig2_train = sum(r.^2) / n; %sum of squared residuals divided by length of model order
                R2_train = ss.rsquare;
                ss = rmfield(ss, 'r');
                ss.R2_train   = R2_train;
                ss.sig2_train = sig2_train;
                stats_linear{i_fold}{i_win, i_k, i_sensor} = ss;
                
                clear ss
            end
        end
    end
    
    %% 11) Characterize model fit on test set
    for i_k = 1:length(ks)
        disp(['ToneDur: ' tonedur_title '0ms - Test Set (Fold: ' ...
            num2str(i_fold) ') - Current k-model: ' num2str(i_k)])
        for i_win = 1:size(windows,1)
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
                
                %Find error of training set on test set
                %Take stats from training set (for current time window,
                %model order, and sensor) and apply them to test set
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
                r = dv - dv_est';
                %compute residuals/difference between estimated and real
                %test set ERF activity for each trial across all tones from 16:32
                
                SS_res = sum(r.^2);%Compute summed (across trials) square of residual
                SS_tot = sum((dv - mean(dv)).^2);%Compute summed square of single trial across-trial average ERF
                R2_test = 1 - SS_res / SS_tot;   %r-square for test set based on estimated data
                
                n = length(r);
                sig2_test = sum(r.^2) / n; %k-value by which we later define the k-prime value
                ss.R2_test   = R2_test;
                ss.sig2_test = sig2_test;
                
                %linear regression on test set
                ss2 = regstats(dv, iv, 'linear',  {'tstat', 'fstat', 'r', 'rsquare', 'adjrsquare'});
                r = ss2.r;
                n = length(r);
                sig2_test_regr = sum(r.^2) / n;
                
                ss.R2_test_regr   = ss2.rsquare;
                ss.sig2_test_regr = sig2_test_regr;
                
                stats_linear{i_fold}{i_win, i_k, i_sensor} = ss; %Add test stats to common stats output var
                
                clear ss
            end
        end
    end
end

%% 12) Save output variables
savefile = [path_save sub '_ExpK_' tonedur_title '.mat'];
save(savefile, 'stats_linear', 'Info_TrialSelection');

end