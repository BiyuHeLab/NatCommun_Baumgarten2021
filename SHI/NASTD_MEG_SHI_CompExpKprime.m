function NASTD_MEG_SHI_CompExpKprime...
    (sub, tonedur_text, ...
    plot_figs, save_figs, ...
    paths_NASTD_MEG)

%Aim: Combine k-vaues across folds and determine winning k-value 
%(i.e., k-prime/k'); Plot Topoplot showing single-subject K-prime values

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

path_inputdata = [paths_NASTD_MEG.Current_outputdata 'ExpK/' tonedur_title '/' sub '/'];
path_outputdata = [paths_NASTD_MEG.Current_outputdata 'ExpK/' tonedur_title '/' sub '/'];
mkdir(path_outputdata)

%% 1) Load k values data for all folds
load([path_inputdata sub '_ExpK_' tonedur_title '.mat']); %var name: stats_linear

%% 2) Define analysis parameters
NumFolds = length(stats_linear);

win_size    = 30;
win_overlap = 0;

samplingFreq = 600;
toneDur_inSecs  = str2num(tonedur_text);
nSamplesPerTone = toneDur_inSecs * samplingFreq;

windows = [1 win_size];
while windows(end,end) < nSamplesPerTone
    windows = [windows; windows(end,:) + (win_size - win_overlap)];
end

if windows(end,end) > nSamplesPerTone
    windows(end,:) = [];
end

windows_inms = (windows / samplingFreq) * 1000;

nSensors = size(stats_linear{1},3);
nWindows = length(windows);
n_k      = size(stats_linear{1},2); %number model order;

k_list = 1:n_k;

%% 3) Determine Kprime as 'winning' k-model (model with minimal sum squared residuals)
%Restructure regstats data and read out relevant stats for model selection from even set
for i_win = 1:nWindows
    for i_sensor = 1:nSensors
        for i_fold = 1:NumFolds
            
            for i_k = 1:n_k
                sumsquaredRes_proxy(i_k) = stats_linear{i_fold}{i_win, i_k, i_sensor}.sig2_test;
                %sum of squared residuals divided by length of model order for each k (per sensor,time win/fold)
                R2_TrainModelonTestData(i_k) = stats_linear{i_fold}{i_win, i_k, i_sensor}.R2_test;
                %R-square (amount of variance in test set explained by the train model) for each k (per sensor,time win/fold)                
            end
            
            %Compute criteria for model selection per set and  place timewin/sensor wise parameter from all folds in common struct            
            [min_sumsqaredRes, index_min_sumsquaredRes] = min( sumsquaredRes_proxy );
            %read out minimum squared residuals divided by length of model
            %order and the respective position in array per sensor/time win/fold
            
            kprime_allFolds{i_win}(i_fold,i_sensor) = k_list(index_min_sumsquaredRes) - 1;
            %k-prime-value (i.e., preferred model order/Number of tones
            %back in sequence history best explaining MEG data;
            %(corrected by -1 because k-prime reflects how many tones back 
            %explains neural activity (i.e., curren tone = 0), whereas 
            %k_list reflects the tone index (i.e., current tone = 1)
            
            kprime_sumsquaredRes_allFolds{i_win}(i_fold,i_sensor) = min_sumsqaredRes;
            %minimum squared residuals divided by length of model order
            %of K-value (i.e., preferred model order/Number of tones back)
            
            sumsquaredRes_allK_allFolds{i_win}{i_sensor}(i_fold, :) = sumsquaredRes_proxy;
            %sum of squared residuals divided by length of model order,
            %stored in common matrix across folds

            R2_TrainModelonTestData_Kprime_allFolds{i_win}(i_fold,i_sensor) = R2_TrainModelonTestData(index_min_sumsquaredRes);
            %R-square (amount of variance in test set explained by the train model) for winning k-model
            %stored in common matrix across folds
            
            R2_TrainModelonTestData_allK_allFolds{i_win}{i_sensor}(i_fold, :) = R2_TrainModelonTestData;
            %R-square (amount of variance in test set explained by the train model) per k-model
            %stored in common matrix across folds
                        
            %Average beta regression weights for each k-model across sets
            for i_k = 1:n_k
                BetaRegWeights_allFolds{i_win, i_k, i_sensor}(i_fold,:) =  ...
                    stats_linear{i_fold}{i_win, i_k, i_sensor}.tstat.beta';
            end
            
            %Clean-up
            clear sumsquaredRes_proxy min_sumsqaredRes index_min_sumsquaredRes R2_TrainModelonTestData
            
        end
    end
end

%% 4) Combine Kprime across folds by averaging values across cross-validation folds
%i.e., average model selection criteria by averaging across folds
for i_win = 1:nWindows
    for i_sensor = 1:nSensors
        
        %Average and compute STD across folds (for each sensor/time win)
        %Kprime
        kprime_AVGacrossFolds{i_win}(i_sensor) = ...
            mean(kprime_allFolds{i_win}(:,i_sensor));
        kprime_STDacrossFolds{i_win}(i_sensor) = ...
            std(kprime_allFolds{i_win}(:,i_sensor));
        %SSR Kprime
        kprime_sumsquaredRes_AVGacrossFolds{i_win}(i_sensor) = ...
            mean(kprime_sumsquaredRes_allFolds{i_win}(:,i_sensor));
        kprime_sumsquaredRes_STDacrossFolds{i_win}(i_sensor) = ...
            std(kprime_sumsquaredRes_allFolds{i_win}(:,i_sensor));
        %SSR all K
        sumsquaredRes_allK_AVGacrossFolds{i_win}(i_sensor, :) = ...
            mean(sumsquaredRes_allK_allFolds{i_win}{i_sensor});
        sumsquaredRes_allK_STDacrossFolds{i_win}(i_sensor, :) = ...
            std(sumsquaredRes_allK_allFolds{i_win}{i_sensor});
        
        %R^2 Kprime
        R2_TrainModelonTestData_Kprime_AVGacrossFolds{i_win}(i_sensor) = ...
            mean(R2_TrainModelonTestData_Kprime_allFolds{i_win}(:,i_sensor));
        R2_TrainModelonTestData_Kprime_STDacrossFolds{i_win}(i_sensor) = ...
            std(R2_TrainModelonTestData_Kprime_allFolds{i_win}(:,i_sensor));
        %R^2 all K
        R2_TrainModelonTestData_allK_AVGacrossFolds{i_win}(i_sensor, :) = ...
            mean(R2_TrainModelonTestData_allK_allFolds{i_win}{i_sensor});
        R2_TrainModelonTestData_allK_STDacrossFolds{i_win}(i_sensor, :) = ...
            std(R2_TrainModelonTestData_allK_allFolds{i_win}{i_sensor});
        
        %Beta Regression Weights
        %all K
        for i_k = 1:n_k
            BetaRegWeights_AVGacrossFolds{i_win, i_k, i_sensor} = ...
                mean(BetaRegWeights_allFolds{i_win, i_k, i_sensor});
            BetaRegWeights_STDacrossFolds{i_win, i_k, i_sensor} = ...
                std(BetaRegWeights_allFolds{i_win, i_k, i_sensor});
        end
        
        %K-prime
        kprime_BetaRegWeight_AVGacrossFolds{i_win}{i_sensor} = ...
            mean(BetaRegWeights_AVGacrossFolds{i_win, ...
            round(kprime_AVGacrossFolds{i_win}(i_sensor) + 1), i_sensor});
        kprime_BetaRegWeight_STDacrossFolds{i_win}{i_sensor} = ...
            std(BetaRegWeights_AVGacrossFolds{i_win, ...
            round(kprime_AVGacrossFolds{i_win}(i_sensor) + 1), i_sensor});
    end
end

%Create summary file containing all relevant vars
Kprime_data = struct;

Kprime_data.perFold = struct;
Kprime_data.perFold.Kprime = kprime_allFolds;
Kprime_data.perFold.SSR_perK = sumsquaredRes_allK_allFolds;
Kprime_data.perFold.SSR_Kprime = kprime_sumsquaredRes_allFolds;
Kprime_data.perFold.BetaRegWeights_perK = BetaRegWeights_allFolds;

Kprime_data.avgFolds = struct;
Kprime_data.avgFolds.Kprime.Avg = kprime_AVGacrossFolds;
Kprime_data.avgFolds.Kprime.STD = kprime_STDacrossFolds;
Kprime_data.avgFolds.SSR_Kprime.Avg = kprime_sumsquaredRes_AVGacrossFolds;
Kprime_data.avgFolds.SSR_Kprime.STD = kprime_sumsquaredRes_STDacrossFolds;
Kprime_data.avgFolds.SSR_perK.Avg = sumsquaredRes_allK_AVGacrossFolds;
Kprime_data.avgFolds.SSR_perK.STD = sumsquaredRes_allK_STDacrossFolds;
Kprime_data.avgFolds.R2_Kprime.Avg = R2_TrainModelonTestData_Kprime_AVGacrossFolds;
Kprime_data.avgFolds.R2_Kprime.STD = R2_TrainModelonTestData_Kprime_STDacrossFolds;
Kprime_data.avgFolds.R2_perK.Avg = R2_TrainModelonTestData_allK_AVGacrossFolds;
Kprime_data.avgFolds.R2_perK.STD = R2_TrainModelonTestData_allK_STDacrossFolds;
Kprime_data.avgFolds.BetaRegWeights_perK.Avg = BetaRegWeights_AVGacrossFolds;
Kprime_data.avgFolds.BetaRegWeights_perK.STD = BetaRegWeights_STDacrossFolds;
Kprime_data.avgFolds.BetaRegWeights_Kprime.Avg = kprime_BetaRegWeight_AVGacrossFolds;
Kprime_data.avgFolds.BetaRegWeights_Kprime.STD = kprime_BetaRegWeight_STDacrossFolds;

%% 5) Save output variables
savefile = [path_outputdata sub '_ExpKprime_' tonedur_title '.mat'];
save(savefile, 'Kprime_data');

%% 6) Plot topoplot showing Kprime values per sensor for each time window
if plot_figs == 1
    %Load in label file
    load([paths_NASTD_MEG.ScriptsDir 'MEG_sensor_setup_272/label272.mat'] ); %file with CTF sensor labels for 272 sensors
    
    %find zlimits for constant plotting across time windows
    dv_absmax = 0;
    dv_absmin = 16;
    
    for i_win = 1:size(windows,1)
        dv_win_absmax = round(max(Kprime_data.avgFolds.Kprime.Avg{i_win}));
        dv_win_absmin = round(min(Kprime_data.avgFolds.Kprime.Avg{i_win}));
        
        if dv_win_absmax > dv_absmax
            dv_absmax = dv_win_absmax;
        end
        if dv_win_absmin < dv_absmin
            dv_absmin = dv_win_absmin;
        end
    end
    
    %Plot topo for current window
    for i_win = 1:size(windows,1)
        samplingFreq = 600;
        %determine start & end points for respective window
        w1 = num2str(1000 * windows(i_win, 1) / samplingFreq , 3 );
        w2 = num2str(1000 * windows(i_win, 2) / samplingFreq , 3 );
        
        model_comp_title = 'Kprime per sensor';
        win_title = ['TW = [' w1 ' ms - ' w2 ' ms]'];
        
        clear dat
        dat.dimord = 'chan_time';
        dat.label  = label;
        dat.time   = 0;
        
        cfg = [];
        cfg.layout    = 'CTF275.lay';
        cfg.comment   = 'no';
        cfg.colorbar  = 'yes';
        cfg.zlim      = [dv_absmin, dv_absmax]; % 'maxabs';
        
        dv = Kprime_data.avgFolds.Kprime.Avg{i_win}'; %read out k values per sensor
        dat.avg = dv;
        
        h = figure;
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        
        ft_topoplotER(cfg, dat);
        title({[sub ' ' model_comp_title], ['TD: ' tonedur_title '0ms; ERF ' win_title]})
        
        if save_figs    == 1
            path_figs = [paths_NASTD_MEG.Current_outputfig 'ExpKprime/' tonedur_title '/' sub '/'];
            mkdir(path_figs)
            
            filename     = [sub '_Kprime_' tonedur_title '_TW' w1 '-' w2 '.png'];
            figfile      = [path_figs filename];
            saveas(gcf, [figfile], 'png'); %save png version
            delete(h);
        end
    end
end

end