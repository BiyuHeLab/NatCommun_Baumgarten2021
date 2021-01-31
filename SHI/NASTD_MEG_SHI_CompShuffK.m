function NASTD_MEG_SHI_CompShuffK...
    (sub, tonedur_text, ...
    NumFolds,nReps, ...
    paths_NASTD_MEG)
%Aim: Compute k-values based on shuffled tone order to create basis for null distribution
%(6-fold cross-validation procedure as linear regression (across trials)
%between sensor-wise, window-wise MEG activity and tone pitch for tones
%32:16)
%For each specified subject and tone duration:
%1) Compute null-hypothesis k-values by randomly permuting tone pitch sequences,
%thus destryoing the temporal dependency of the tone pitches
%Do tone order shuffling identically for each tone duration
%condition, to allow for valid comparison across tone duration conditions

%% 0) Specify vars, paths, and setup fieldtrip
addpath(genpath(paths_NASTD_MEG.ScriptsDir));
NASTD_MEG_SubInfo

%Data input
loadfile_MEG = [si.path_ica  'data_clean.mat'];%path to single-subject MEG data
loadfile_behav = si.path_behav;%path to single-subject behavioral/stimulus data
 
%Data output for Info file
path_save = [paths_NASTD_MEG.Current_outputdata 'ShuffK/'];
mkdir(path_save)

%% 1) load preprocessed MEG data
tic
disp('loading...')
load(loadfile_MEG);
load(loadfile_behav);
disp(['done loading in ' num2str(toc) ' sec'])

data_raw = data; %backup data because subsequently it is reduced to specific tone dur

for i_tonedur = 1:length(tonedur_text)
    
    %Tone duration condition for data load-in
    if strcmp(tonedur_text{i_tonedur},'0.15')
        tonedur_title = '0.15sTD';
        tonedur_field = 'short';
    elseif  strcmp(tonedur_text{i_tonedur},'0.3')
        tonedur_title = '0.3sTD';
        tonedur_field = 'med';
    elseif strcmp(tonedur_text{i_tonedur},'0.6')
        tonedur_title = '0.6sTD';
        tonedur_field = 'long';
    end
        
    %Data output
    path_save_ShuffK = [paths_NASTD_MEG.Current_outputdata 'ShuffK/' tonedur_title '/' sub '/'];
    mkdir(path_save_ShuffK)
    
    %% 2) Bring data into required format
    data = data_raw;
    
    %Extract trials from clean MEG data
    data_MEG.(tonedur_field) = extract_trials(data_clean, sub, str2double(tonedur_text{i_tonedur}));
    
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
    dat_MEG_fields  = fieldnames(data_MEG.(tonedur_field));
    for i = 1:length(dat_MEG_fields)
        if eval(['length(data_MEG.(tonedur_field).' dat_MEG_fields{i} ') ==  ' numoftrials])
            eval(['data_MEG.(tonedur_field).' dat_MEG_fields{i} ' = data_MEG.(tonedur_field).' dat_MEG_fields{i} '(si.good_trials);']);
        end
    end
    
    %Determine and label individual sequence ID
    %9 unique sequences (3 beta * 3 p*34)* 3 tone dur = 27 trials/blocks, 12
    %blocks total = 324 trials total
    %Create seqID proxy field
    data.stim.seqID = (1:length(data.stim.series_f));
    data.stim.seqID(1:end) = NaN;
    seqID_perUniqueSum = [];
    
    %Read out all trials of first block and determine individual sequence ID
    for i_trial = 1:length(data.stim.series_f)/12 %trials per block (27)
        All_trials_block1(i_trial,:) = data.stim.series_f{i_trial}(1:33); %only first 33 tone, since 34th varies
    end
    seqID_perUniqueSum = unique(sum(All_trials_block1')); %use sum across tone frequency values as metric determinig SeqID
    %should result in 9 unique sequences (independent of final tone pitch/p34)
    %Identical across tone duration conditions
    
    %Individually label all trials based on their across tone freq sum
    for i_trial = 1:length(data.stim.series_f) %loop across all trials
        for i_seqID = 1:length(seqID_perUniqueSum) %loop across uSeqs
            if sum(data.stim.series_f{i_trial}(1:33)) == seqID_perUniqueSum(i_seqID)
                data.stim.seqID(i_trial) = i_seqID;
            end
        end
    end
    
    %Filter behavioral data & MEG data for trials corresponding to a particular tone length
    numoftrials = num2str(length(data.trialNum));
    filter_tonedur = data.stim.toneDur == str2double(tonedur_text{i_tonedur});
    
    %Filter data
    dat_fields  = fieldnames(data);
    for i = 1:length(dat_fields)
        if eval(['length(data.' dat_fields{i} ') == ' numoftrials])
            eval(['data.' dat_fields{i} ' = data.' dat_fields{i} '(filter_tonedur);']);
        end
    end
    
    %filter data.stim
    stim_fields = fieldnames(data.stim);
    for i = 1:length(stim_fields)
        if eval(['length(data.stim.' stim_fields{i} ') == ' numoftrials])
            eval(['data.stim.' stim_fields{i} ' = data.stim.' stim_fields{i} '(filter_tonedur);']);
        end
    end
    
    %filter data_MEG
    MEG_fields = fieldnames(data_MEG.(tonedur_field));
    for i = 1:length(MEG_fields)
        if eval(['length(data_MEG.(tonedur_field).' MEG_fields{i} ') == ' numoftrials])
            eval(['data_MEG.(tonedur_field).' MEG_fields{i} ' = data_MEG.(tonedur_field).' MEG_fields{i} '(filter_tonedur);']);
        end
    end
    
    %Tone Duration Filtered Sanity check: count number of seqID trials 
    %(should be 12 per seqID per tone duration condition)
    for i_seqID = 1:length(seqID_perUniqueSum)
        if length(find(data.stim.seqID==i_seqID)) ~= 12
            warning0 = ['CAVE: IMBALANCE IN NUMBER OF REPS FOR uSeq:' num2str(i_seqID) ' IN TONEDUR CONDITION ' num2str(i_tonedur)];
            disp(warning0);
        end
    end
    
    NaNtrial = [];
    %Check MEG data for NaNs
    for i_trial = 1:length(data_MEG.(tonedur_field).trial)
        NaNfind_dataMEG{i_trial} = find(isnan(data_MEG.(tonedur_field).trial{i_trial}));
        if ~isempty(NaNfind_dataMEG{i_trial})
            NaNtrial(1,:) = i_trial;
        end
    end
    
    if ~isempty(NaNtrial)
        newTrialSelection = 1:length(data_MEG.(tonedur_field).trial);
        newTrialSelection(NaNtrial) = [];
        
        sampleinfo = data_MEG.(tonedur_field).sampleinfo;
        sampleinfo(NaNtrial) = [];
        %If NaNs present, delete trials with NaNs
        cfg = [];
        cfg.trials =  newTrialSelection;
        data_MEG.(tonedur_field) = ft_selectdata(cfg, data_MEG.(tonedur_field));
        data_MEG.(tonedur_field).sampleinfo = sampleinfo;
    end
    
    %% 3) Define folds, allocate trials to folds
    %Importantly: Balanced split per unique tone seqence (12 repetitions for
    %each of 9 unique tone sequences (1:32) for each of 3 tone duration conditons)
    %Final tones (p34) are not taken into account here, since history tracking
    %analysis focusses only on tone 1:32
    
    %Determine uSeq-wise trial selection
    Number_uSeq = length(unique(data.stim.seqID));
    
    %Create empty proxies
    Number_reps_peruSeq = NaN(1,Number_uSeq);
    index_reps_peruSeq = NaN(Number_uSeq,12);
    
    index_reps_peruSeq_noNaN = [];
    
    index_train = [];
    index_test = [];
    
    %Find trial indices corresponding to reps in each uSeq
    for i_seqID = 1:Number_uSeq
        %find index/position of all uSeq reps in tonedur specific data struct
        index_reps_peruSeq(i_seqID,1:length(find(data.stim.seqID == i_seqID))) = find(data.stim.seqID == i_seqID);
        
        %Check if any trial (determined by data.stim) is not present in MEG data
        for i_entry = 1:length(index_reps_peruSeq(i_seqID,1:length(find(data.stim.seqID == i_seqID))))
            if index_reps_peruSeq(i_seqID,i_entry) > length(data_MEG.(tonedur_field).trial)
                index_reps_peruSeq(i_seqID,i_entry) = NaN;
            end
        end
        
        %get rid of NaN entries (requires cell due to potential different
        %vector length between uSeqs)
        counter = 1;
        for i_rep = 1:length(index_reps_peruSeq)
            if ~isnan(index_reps_peruSeq(i_seqID,i_rep))
                index_reps_peruSeq_noNaN{i_seqID}(counter)  = index_reps_peruSeq(i_seqID,i_rep);
                counter = counter+1;
            end
        end
    end
    
    %Allocate reps to groups/folds (optimal if identical number of reps per
    %group, but not necessary)
    
    %Update Number of folds based on reps per uSeq (to guarantee that no
    %fold is empty
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
    
    %Shuffle reps
    for i_seqID = 1:Number_uSeq
        index_reps_peruSeq_noNaN_shuffled{i_seqID} = Shuffle(index_reps_peruSeq_noNaN{i_seqID});
        
        %count used repeptitions per unique tone sequence(1:33)
        Info_TrialSelection.(tonedur_field).Number_Usedreps_peruSeq(1,i_seqID) = ...
            length(index_reps_peruSeq_noNaN{i_seqID});
        Info_TrialSelection.(tonedur_field).TrialIndexUsedreps_ToneDurDataSet_peruSeq{i_seqID} = ...
            index_reps_peruSeq_noNaN{i_seqID};
        Info_TrialSelection.(tonedur_field).TrialIndexUsedreps_WholeDataSet_peruSeq{i_seqID} = ...
            data.trialNum(index_reps_peruSeq_noNaN{i_seqID});
        
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
            %For uneven number of reps, check that allocated reps are present
            %in data set, otherwise delete unpresent one (leads to different
            %number of uSeq-reps between folds - fewer reps placed in final fold)
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
        Info_TrialSelection.(tonedur_field).Number_Usedreps_perFold{i_fold} = ...
            length(index_reps_peruSeq_perfold{i_fold});
        Info_TrialSelection.(tonedur_field).TrialIndexUsedreps_ToneDurDataSet_perFold{i_fold} = ...
            index_reps_peruSeq_perfold{i_fold};
        Info_TrialSelection.(tonedur_field).TrialIndexUsedreps_WholeDataSet_perFold{i_fold} = ...
            data.trialNum(index_reps_peruSeq_perfold{i_fold});
        
        %Store which reps/trials per uSeq were allocated to which fold
        for i_seqID = 1:Number_uSeq
            Info_TrialSelection.(tonedur_field).TrialIndexUsedreps_ToneDurDataSet_uSeqRepsperFold{i_fold}{i_seqID} = ...
                intersect(index_reps_peruSeq_perfold{i_fold},index_reps_peruSeq_noNaN{i_seqID});
            Info_TrialSelection.(tonedur_field).TrialIndexUsedreps_WholeDataSet_uSeqRepsperFold{i_fold}{i_seqID} = ...
                data.trialNum(intersect(index_reps_peruSeq_perfold{i_fold},index_reps_peruSeq_noNaN{i_seqID}));
        end
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
        f_train = [1:length(data_MEG.(tonedur_field).trial)]*0; %create empty filter-proxies
        f_test = [1:length(data_MEG.(tonedur_field).trial)]*0;
        
        f_train(index_train) = 1; %fill set-specific chosen trial indices
        f_test(index_test) = 1;
        
        f_train = logical(f_train); %convert arrays to logical
        f_test = logical(f_test);
                
        %% 5) Define training set
        %Copy field struct and input from selected train trials
        data_MEG_train = data_MEG.(tonedur_field);
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
        
        clear data_train
        
        %% 6) Define test set
        %Copy field struct and input from selected train trials
        data_MEG_test = data_MEG.(tonedur_field);
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
        
        clear data_test
        
        %% 7) Analyse time window defined MEG data
        %Define analysis parameters
        nTrials_train = length(data_MEG_train.trial);
        nTrials_test  = length(data_MEG_test.trial);
        nSensors = size( data_MEG_train.trial{1}, 1 );
        
        samplingFreq              = data_MEG_train.fsample;
        toneDur_inSecs  = tonedur_text{i_tonedur};
        nSamplesPerTone = str2double(toneDur_inSecs) * samplingFreq;
        
        nSamplesPerSeries = length(data_MEG_train.trial{1}(1,:));
        
        %Define time windows used for analysis
        win_size    = 30; 
        win_overlap = 0; 
        
        %Define starting point (i.e., first tone)
        t_start_ind(1) = find(data_MEG_train.time{1} == 0); %Start of first tone defined as t = 0
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
        
        if windows(end,end) > nSamplesPerTone %TJB: windows within single tone
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
            tone_pitch_train(i_trial, :) = log( data_MEG_train.stim.series_f{i_trial} );
            %output: trial*tone matrix containing tone pitches (freq) of all tones per sequence
        end
        
        for i_trial = 1:nTrials_test
            tone_pitch_test(i_trial, :)  = log( data_MEG_test.stim.series_f{i_trial} );
        end
        
        tone_pitch_train = tone_pitch_train(:, [1:end-2]);
        %take original matrix, delete last two column/tone, because we only use first 32 tone pitches
        tone_pitch_test = tone_pitch_test(:, [1:end-2]);
        
        %% 10) Shuffle temporal order of tone pitches (tone 1:32) in training set
        %Shuffle details:
        %1) Shuffle tone pitch order anew for every run
        %2) Across runs and tone durations, shuffle tone pitch order
        %the same way for all identical uSeqs (i.e., not for every
        %trial, but for every uSeq defined by 1:33 tone)
        %3) Shuffle tone pitch order anew for every fold
        
        tone_pitch_train_orig = tone_pitch_train; %store original tone pitch order
        
        for i_rep = 1:nReps %loop for shuffle repetitions
            disp(['--- Current Rep# = ' num2str(i_rep) ' ---'])
            clear tone_pitch_train %delete var for each Rep-run
            if i_tonedur == 1 %for tonedur = 150 ms, shuffle and record shuffle order for each uSeq/trial
                for i_uSeq = unique(data_MEG_train.stim.seqID) %for each uSeq in train set
                    uSeq_trials = find(data_MEG_train.stim.seqID == i_uSeq); %find indices of all trials corresponding to current uSeq
                    t(1,:) = tone_pitch_train_orig(uSeq_trials(1),1:32); %save original tone order for current uSeq (taken from first trial corresponding to current uSeq)
                    t(2,:) = 1:32; %Create label line to track shuffling order
                    
                    %%Critical Step: Shuffling of original tone order for current uSeq
                    t = Shuffle(t,1); %shuffle tone order for current uSeq - Shuffle.m function from Psychtoolbox3 - shuffle on 1st dim (rows)
                    for i_trialuSeq = uSeq_trials %apply this shuffling also for the other trials with current uSeq
                        tone_pitch_train(i_trialuSeq,:) = t(1,:);
                        Info_TrialSelection.(tonedur_field).Shuffleorder_allTrials_perFold{i_fold}{i_rep}(i_trialuSeq,:) = t(2,:);  %save shuffle order for each trial
                    end
                    Info_TrialSelection.(tonedur_field).Shuffleorder_peruSeq_perFold{i_fold}{i_rep}(i_uSeq,:) = t(2,:);  %save shuffle order for each uSeq
                end
                
            else %for tonedur = 300  or 600 ms, take shuffling order from 150ms condition and order tones of corresponind uSeq accordingly
                %Load trialinfo output-file containing shuffle order
                load([path_save sub '_ShuffK__Shuffleorder_peruSeq_perFold.mat']);                
                for i_uSeq = unique(data_MEG_train.stim.seqID) %for each uSeq in train set
                    uSeq_trials = find(data_MEG_train.stim.seqID == i_uSeq); %find indices of all trials corresponding to current uSeq
                    for i_trialuSeq = uSeq_trials %apply prior determined shuffling to other trials with current uSeq
                        for i_tone = 1:32
                            tone_pitch_train(i_trialuSeq,i_tone) = ...
                                tone_pitch_train_orig(i_trialuSeq,...
                                Shuffleorder_peruSeq_perFold{i_fold}{i_rep}(i_uSeq,i_tone));
                            %replace current tone of current uSeq-trial
                            %with the uSeq-specific shuffled tone
                            %determined by the first tone dur condition
                            Info_TrialSelection.(tonedur_field).Shuffleorder_peruSeq_perFold{i_fold}{i_rep}(i_uSeq,i_tone) = ...
                                Shuffleorder_peruSeq_perFold{i_fold}{i_rep}(i_uSeq,i_tone);%save shuffle order for each uSeq
                        end
                    end
                end
            end
            
            %% 11) Compute k values based on traning set
            %compute k as linear regression (across trials) between sensor-wise, window-wise ERF and
            %tone pitch for tones 32:16 (17 tones total)
            ks = 1:16; %1:16 regression terms/k-values/models - all potential k vals
            
            for i_k = 1:length(ks) %index/current k val (i.e., amount of k models)
                disp(['ToneDur: ' tonedur_title '0ms - Training Set (Fold: ' num2str(i_fold) ') - Current k-model: ' num2str(i_k)])
                
                for i_win = 1:size(windows,1)
                    for i_sensor = 1:nSensors
                        
                        clear dv iv_cell
                        dv = [];
                        for toneIndex = 16:32
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
                        stats_linear{i_win, i_k, i_sensor} = ss;
                        clear ss
                    end
                end
            end
                        
            %% 12) Characterize model fit on test set
            %compute k as linear regression (across trials) between sensor-wise, window-wise ERF and
            %tone pitch for tones 32:16 (17 tones total), similar to above
            for i_k = 1:length(ks)
                disp(['ToneDur: ' tonedur_title '0ms - Test Set (Fold: ' num2str(i_fold) ') - Current k-model: ' num2str(i_k)])
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
                        %Take stats from training set (for current time window, model order, and sensor)
                        %and apply them to test set
                        ss = stats_linear{i_win, i_k, i_sensor};
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
                        SumSquareRes_Shuff{i_fold}{i_win, i_k, i_sensor}(i_rep) = sig2_test;
                        %output: matrix with sum of squared residuals divided by amount of tested models
                        %per time window, model order, sensor, repetition
                    end
                end
            end
        end
        %Cleanup before each fold
        clear tone_pitch_train tone_pitch_test erf_win_train erf_win_test dv_est
    end
    
    %% 13) Save output variables
    %Save Shuffle tone order to be used by next tone duration conditions
    if i_tonedur == 1 %for tonedur = 150 ms, save shuffle order
        Shuffleorder_peruSeq_perFold = Info_TrialSelection.(tonedur_field).Shuffleorder_peruSeq_perFold;
        savefile = [path_save sub '_ShuffK__Shuffleorder_peruSeq_perFold.mat'];
        save(savefile,'Shuffleorder_peruSeq_perFold');
    end
    
    %Save output var with shuffled SumSquareRes
    savefile = [path_save_ShuffK sub '_ShuffK_' tonedur_title '.mat'];
    save(savefile, 'SumSquareRes_Shuff', 'Info_TrialSelection');
    
    %Cleanup
    clearvars -except ...
        sub si nReps  NumFolds tonedur_text ...  %Input info
        data_raw data_clean ...   %behav and MEG data
        i_tonedur i_fold...   %loop info
        loadfile* path* ... %path info
        Info_TrialSelection ... %shuffle order info
        paths_NASTD_MEG
end

end