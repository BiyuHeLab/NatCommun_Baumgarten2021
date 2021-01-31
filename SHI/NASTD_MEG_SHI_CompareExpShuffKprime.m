function NASTD_MEG_SHI_CompareExpShuffKprime...
    (sub, tonedur_text, ...
    plot_figs, pval_plotting, save_figs, ...
    paths_NASTD_MEG)
%Aim: Combine k-prime values from experimental vs. shuffled null distribution
%1) Combine shuffled k-values across sets by averaging each rep
%across folds to receive null distribution of 100 K-prime values per time
%window/sensor
%2) Compare to experimental k-prime distribution and figure out position of 
%exp. relative to shuffled data
%3) Save combined output
%4) Plot topoplot of sign. Kprime values (exp. vs. shuffle)

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

path_input_Shuffk = [paths_NASTD_MEG.Current_outputdata 'ShuffK/' tonedur_title '/' sub '/'];%path were shuffled K-value output is
path_input_Expk = [paths_NASTD_MEG.Current_outputdata 'ExpK/' tonedur_title '/' sub '/']; %path were exp K-value output is
path_outputdata = [paths_NASTD_MEG.Current_outputdata 'ExpK/' tonedur_title '/' sub '/'];
mkdir(path_outputdata)

%% 1)Load in shuffled Kvalues
load([path_input_Shuffk sub '_ShuffK_' tonedur_title '.mat']);

%Define analysis parameters
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

nSensors = size(SumSquareRes_Shuff{1},3);
nWindows = length(windows);
n_k      = size(SumSquareRes_Shuff{1},2); %number model order;
nReps = length(SumSquareRes_Shuff{1}{1,1,1});

k_list = 1:n_k;

%Process shuffled sets
for i_fold = 1:length(SumSquareRes_Shuff)
    for i_rep = 1:nReps
        for i_win = 1:nWindows
            for i_sensor = 1:nSensors
                for i_k = 1:n_k
                    s(i_k) = SumSquareRes_Shuff{i_fold}{i_win, i_k, i_sensor}(i_rep);
                    %read out shuffled sum of squared residuals divided by length of model order, per model order
                end
                
                [min_s, min_s_ind] = min( s );
                %read out minimum squared residuals divided by length of model order and the respective position in array
                
                Kprime_Shuff.Kprime_perWinSensFoldRep{i_win}{i_sensor}(i_fold, i_rep) = k_list(min_s_ind) - 1;
                %Kprime (i.e., preferred model order/Number of tones back best explaining MEG data due to minimal summed squared residuals)
                %corrected by -1 to correct for lin reg offset term)
                
                Kprime_Shuff.MinSumSquaredRes_perWinSensFoldRep{i_win}{i_sensor}(i_fold, i_rep) = (min_s);
                %minimal summed squared residuals defining Kprime
                
                clear s min_s min_s_ind %cleanup
            end
        end
    end
end

clear SumSquareRes_Shuff %cleanup

for i_win = 1:nWindows
    for i_sensor = 1:nSensors
        Kprime_Shuff.Kprime_perWinSensFold_AvgRep{i_win}(:,i_sensor) = ...
            mean(Kprime_Shuff.Kprime_perWinSensFoldRep{i_win}{i_sensor},2);
        %average across reps within fold, to later compare avg shuffled
        %Kprime across folds
    end
end


%% 2) Combine all folds for shuffled null distribution by averaging across folds
%Output: distribution of Kprime for 100 reps (averaged across folds) for
%each time window/sensor
for i_win = 1:nWindows
    for i_sensor = 1:nSensors
        Kprime_Shuff.Kprime_perWinSens.NullDist_AvgFolds{i_win}(i_sensor,:) = ...
            mean(Kprime_Shuff.Kprime_perWinSensFoldRep{i_win}{i_sensor});
        %For each rep, average Kprime across folds, distribution of 100 reps
        
        Kprime_Shuff.Kprime_perWinSens.AvgFoldsReps{i_win}(i_sensor,:) = ...
            mean(mean(Kprime_Shuff.Kprime_perWinSensFoldRep{i_win}{i_sensor})); %Average Kprime across folds, then avg across reps
        Kprime_Shuff.Kprime_perWinSens.StdReps{i_win}(i_sensor,:) = ...
            std(mean(Kprime_Shuff.Kprime_perWinSensFoldRep{i_win}{i_sensor})); %Std Kprime across reps
        
        Kprime_Shuff.SumSquareRes_perWinSens.AvgFoldsReps{i_win}(i_sensor,:) = ...
            mean(mean(Kprime_Shuff.MinSumSquaredRes_perWinSensFoldRep{i_win}{i_sensor}));
        Kprime_Shuff.SumSquareRes_perWinSens.StdReps{i_win}(i_sensor,:) = ...
            std(mean(Kprime_Shuff.MinSumSquaredRes_perWinSensFoldRep{i_win}{i_sensor}));
    end
end

%% 3) Load in combined experimental k-prime values
load([path_input_Expk sub '_ExpKprime_' tonedur_title '.mat']);

%Computes the amount of shuffled k-prime values that are
%equal or bigger to the experimental/original k-prime value and
%divides this by the number of reps for each time window and sensr, in each repetition
for i_win = 1:nWindows
    for i_sensor = 1:nSensors
        Kprime_data.avgFolds.Kprime.pval_ExpvsShuff{i_win}(i_sensor) = ...
            sum( Kprime_Shuff.Kprime_perWinSens.NullDist_AvgFolds{i_win}(i_sensor, :) >= ...
            Kprime_data.avgFolds.Kprime.Avg{i_win}(i_sensor) ) ...
            / length(Kprime_Shuff.Kprime_perWinSens.NullDist_AvgFolds{i_win}(i_sensor, :));
        %Add binary significance coding for plotting
        if Kprime_data.avgFolds.Kprime.pval_ExpvsShuff{i_win}(i_sensor) < 0.05
            Kprime_data.avgFolds.Kprime.sign{i_win}(i_sensor) = 1;
        else
            Kprime_data.avgFolds.Kprime.sign{i_win}(i_sensor) = NaN;
        end
    end
end

%% 4) Create summary file with necessary EXP and SHUFF Kprime info (only across-fold averages)
Kprime = struct;
Kprime.Exp.avgFolds = Kprime_data.avgFolds;
Kprime.Shuff.avgFolds.NullDist = Kprime_Shuff.Kprime_perWinSens.NullDist_AvgFolds;
Kprime.Shuff.avgFoldsReps.Kprime.Avg = Kprime_Shuff.Kprime_perWinSens.AvgFoldsReps;
Kprime.Shuff.avgFoldsReps.Kprime.STD = Kprime_Shuff.Kprime_perWinSens.StdReps;
Kprime.Shuff.avgFoldsReps.minSSR.Avg = Kprime_Shuff.SumSquareRes_perWinSens.AvgFoldsReps;
Kprime.Shuff.avgFoldsReps.minSSR.STD = Kprime_Shuff.SumSquareRes_perWinSens.StdReps;

%% 5) Save output variables
%Combined output file with
%1) combined experimental k-prime values,
%2) combined shuffled k-prime values,
savefile = [path_outputdata sub '_ExpvsShuffKprime_' tonedur_title '.mat' ];
save(savefile, 'Kprime');

%% 6) Plot topoplot showing significant (exp. vs. shuff_ Kprime values per sensor for each time window
if plot_figs == 1
    %Load in label file
    load([paths_NASTD_MEG.ScriptsDir 'MEG_sensor_setup_272/label272.mat'] ); %file with CTF sensor labels for 272 sensors
    
    %Find zlimits for constant plotting across time windows
    dv_absmax = 0;
    dv_absmin = 16;
    
    for i_win = 1:size(windows,1)
        dv_win_absmax = round(max(Kprime.Exp.avgFolds.Kprime.Avg{i_win})); 
        dv_win_absmin = round(min(Kprime.Exp.avgFolds.Kprime.Avg{i_win})); 
        
        if dv_win_absmax > dv_absmax
            dv_absmax = dv_win_absmax;
        end
        if dv_win_absmin < dv_absmin
            dv_absmin = dv_win_absmin;
        end
    end
        
    %Plot topo for current window
    for i_win = 1:size(windows,1)
        
        %determine start & end points for respective window
        w1 = num2str(1000 * windows(i_win, 1) / samplingFreq , 3 );
        w2 = num2str(1000 * windows(i_win, 2) / samplingFreq , 3 );
        
        model_comp_title = 'Exp Kprime per sensor';
        title_text = ['p < ' num2str(pval_plotting) ' (uncorrected)'];
        win_title = ['TW = [' w1 ' ms - ' w2 ' ms]'];
        
        clear dat
        dat.dimord = 'chan_time';
        dat.label  = label;
        dat.time   = 0;
        
        cfg = [];
        cfg.layout    = 'CTF275.lay';
        cfg.comment   = 'no';
        cfg.colorbar  = 'yes';
        cfg.zlim      = [dv_absmin, dv_absmax];
        
        cfg.highlight = 'on';
        cfg.highlightchannel = label(Kprime.Exp.avgFolds.Kprime.pval_ExpvsShuff{i_win} < pval_plotting);
        cfg.highlightsymbol  = '.';
        cfg.highlightcolor   = [0.9 0.9 0.9];
        cfg.highlightsize    = 30;
        
        dv = Kprime.Exp.avgFolds.Kprime.Avg{i_win}'; %read out k values per sensor
        dat.avg = dv;
        
        h = figure;
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        
        ft_topoplotER(cfg, dat);
        title({[sub ' ' model_comp_title], ['TD: ' tonedur_title '; ERF ' win_title '; ' , title_text]})
        
        if save_figs    == 1
            path_figs = [paths_NASTD_MEG.Current_outputfig 'ExpKprime/' tonedur_title '/' sub '/'];
            mkdir(path_figs)
            
            filename     = [sub '_SignKprime_' tonedur_title '0_TW' w1 '-' w2 '.png'];
            figfile      = [path_figs filename];
            saveas(gcf, [figfile], 'png'); %save png version
            delete(h);
        end
    end    
end

end