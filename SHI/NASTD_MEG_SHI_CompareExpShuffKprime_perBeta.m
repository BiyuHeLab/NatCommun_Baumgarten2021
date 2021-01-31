function NASTD_MEG_SHI_CompareExpShuffKprime_perBeta...
    (sub, betalevel_input, tonedur_text,...
    NumFolds, nReps, ...
    paths_NASTD_MEG)

%Aim: Combine EXP and SHUFF Kprime results to create Ssub basis for 
%statistical Kprime evaluation
%Implementation:
%1) Combine null distribution k-values across sets by averaging each rep
%across folds to receive null distribution of 100 K-prime values per time
%window/sensor
%2) Add experimental distribution and figure out position of exp. relative
%to shuffled data
%3) Save combined output

%% 0) specify paths and set up vars
addpath(genpath(paths_NASTD_MEG.ScriptsDir));
NASTD_MEG_subjectinfo

%Data input/output
path_load_EXPk = [paths_NASTD_MEG.Current_outputdata 'ExpK/SeqBeta_' betalevel_input '/'];
path_load_SHUFFk = [paths_NASTD_MEG.Current_outputdata 'ShuffK/SeqBeta_' betalevel_input '/'];
path_save = [paths_NASTD_MEG.Current_outputdata 'ExpK/SeqBeta_' betalevel_input '/'];
mkdir(path_save); %make new directory for output files

%% 1)Load in Kprime shuffled null distribution (all folds in 1 file)
%Load set (with identical tone order shuffling)
load([path_load_SHUFFk sub '_ShuffK_SeqBeta' betalevel_input '.mat']);

%Define analysis parameters
nSensors = size(SumSquareRes_Shuff{1},3);
nWindows = size(SumSquareRes_Shuff{1},1);
n_k      = size(SumSquareRes_Shuff{1},2); %number model order;
nReps = length(SumSquareRes_Shuff{1}{1,1,1});

k_list = 1:n_k;

%Process k-fold sets
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
                
                Kprime_Shuff.Kprime_perWinSensFoldRep{i_win}{i_sensor}(i_fold, i_rep) ...
                    = k_list(min_s_ind) - 1;                
                %Kprime (i.e., preferred model order/Number of tones back 
                %best explaining MEG data due to minimal summed squared residuals)

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
%Output: distribution of Kprime for 100 reps (averaged across folds) for each time window/sensor
for i_win = 1:nWindows
    for i_sensor = 1:nSensors    
        Kprime_Shuff.Kprime_perWinSens.NullDist_AvgFolds{i_win}(i_sensor,:) = ...
            mean(Kprime_Shuff.Kprime_perWinSensFoldRep{i_win}{i_sensor}); 
        %For each rep, average Kprime across folds = distribution of 100 reps

        Kprime_Shuff.Kprime_perWinSens.AvgFoldsReps{i_win}(i_sensor,:) = ...
            mean(mean(Kprime_Shuff.Kprime_perWinSensFoldRep{i_win}{i_sensor})); 
        %Average Kprime across folds, then avg across reps    
        Kprime_Shuff.Kprime_perWinSens.StdReps{i_win}(i_sensor,:) = ...
            std(mean(Kprime_Shuff.Kprime_perWinSensFoldRep{i_win}{i_sensor})); 
        %Std Kprime across reps  
        
        Kprime_Shuff.SumSquareRes_perWinSens.AvgFoldsReps{i_win}(i_sensor,:) = ...
            mean(mean(Kprime_Shuff.MinSumSquaredRes_perWinSensFoldRep{i_win}{i_sensor})); 
        Kprime_Shuff.SumSquareRes_perWinSens.StdReps{i_win}(i_sensor,:) = ...
            std(mean(Kprime_Shuff.MinSumSquaredRes_perWinSensFoldRep{i_win}{i_sensor}));         
    end
end

%% 3) Load in combined unshuffled k-prime values and compute number of shuffled k-prime higher than original
%Load in combined unshuffled k-prime values
load([path_load_EXPk  sub '_ExpKprime_SeqBeta' betalevel_input '.mat']);

%for each time window and sensr, in each repetition
for i_win = 1:nWindows    
    for i_sensor = 1:nSensors

        Kprime_data.avgFolds.Kprime.pval_ExpvsShuff{i_win}(i_sensor) = ...
            sum( Kprime_Shuff.Kprime_perWinSens.NullDist_AvgFolds{i_win}(i_sensor, :) >= ...
            Kprime_data.avgFolds.Kprime.Avg{i_win}(i_sensor) ) ...
            / length(Kprime_Shuff.Kprime_perWinSens.NullDist_AvgFolds{i_win}(i_sensor, :));
        %i.e., computes the amount of shuffled k-prime values that are
        %equal or bigger to the experimental/original k-prime value and
        %divides this by the number of reps
        
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
%Combined output file with 1) combined original k-prime values, combined
%2) shuffled k-prime values, 3) amount of shuffled_
savefile = [path_save sub '_ExpvsShuffKprime_SeqBeta' betalevel_input '.mat' ];
save(savefile, 'Kprime');
           
end