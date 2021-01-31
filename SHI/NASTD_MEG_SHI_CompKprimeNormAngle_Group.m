function NASTD_MEG_SHI_CompKprimeNormAngle_Group...
    (param, ...
    paths_NASTD_MEG)
%Aim: Compute Vector Norm, Angle to Duration Line, Angle to Inforation line for each sensor

%% 0) Specify vars, paths, and setup fieldtrip
addpath(genpath(paths_NASTD_MEG.ScriptsDir));

%Data input
path_inputdata  = [paths_NASTD_MEG.Current_analysis 'SummaryStruct/'];

%Data output
path_outputdata = [paths_NASTD_MEG.Current_outputdata 'ExpK/TDcomparison/'];
mkdir(path_outputdata)

%% 1) Load input data
%Input data: SSub & GAvg EXP & SSub SHUFF k-prime vals from summary struct
load([path_inputdata 'KprimeSummaryStruct_SSubsGAvg.mat']);

%% 2. Restructure input data
%Experimental k-prime values (Group-level)
kExpData = Kprime_GAvg.Exp; %across-subject-averaged experimental k-prime values;

% Shuffled k-prime values (null distribution)
%Create new shuffled data struct with k-prime vals for all subjects
for i_tonedur = 1:size(Kprime_AllSub.Shuff,2)
    for i_win = 1:size(Kprime_AllSub.Shuff{1},2) %3 Time windows
        for i_sub = 1:size(Kprime_AllSub.Shuff,1)
            kKeepShuffles_allSubs{i_win,i_tonedur}(i_sub,:,:) = Kprime_AllSub.Shuff{i_sub,i_tonedur}{i_win};
            %Appends k-prime vals per timewin/tonedur for all subs (subs/sensors/reps)
            %Output: X k-prime values for each sensor per subject/timewin/tonedur
        end
        
        %Create group-level null distribution
        %A 1000 times, randomly select a repetition from each subject 
        %with replacement, new selection for each subject)
        for i_samp = 1:param.KprimeComparison.NumSamples4ShuffGroupDist
            for i_sub = 1:size(Kprime_AllSub.Shuff,1)
                concatSubsRand(i_sub, :, i_samp) = squeeze(kKeepShuffles_allSubs{i_win,i_tonedur}...
                    (i_sub, :, randi([1 size(kKeepShuffles_allSubs{1,1},3)])));
            end
        end
        
        %Then average across subs to get distribution of nSamples shuffles
        kNullDistrib{i_win,i_tonedur} = squeeze(mean(concatSubsRand));

    end
end

%Create summary struct for output
Kprime4AllTD = struct;
%Restructure experimental and null distributions to coordinate format
for i_win = 1:size(Kprime_AllSub.Shuff{1},2) %3 Time windows
    
    avg_150ms{1,i_win} = kExpData{1}{i_win};
    avg_300ms{1,i_win} = kExpData{2}{i_win};
    avg_600ms{1,i_win} = kExpData{3}{i_win};
    
    null_150ms{1, i_win} = kNullDistrib{i_win,1};
    null_300ms{1, i_win} = kNullDistrib{i_win,2};
    null_600ms{1, i_win} = kNullDistrib{i_win,3};
    
    x_exp = avg_150ms{1,i_win};
    y_exp = avg_300ms{1,i_win};
    z_exp = avg_600ms{1,i_win};
    
    x_null = null_150ms{1,i_win};
    y_null = null_300ms{1,i_win};
    z_null = null_600ms{1,i_win};
    
    %Store kprime vals for each sensor and ToneDur/Axis in summary struct
    Kprime4AllTD.Coord.perSens.x_exp{i_win} = x_exp;
    Kprime4AllTD.Coord.perSens.y_exp{i_win} = y_exp;
    Kprime4AllTD.Coord.perSens.z_exp{i_win} = z_exp;
    Kprime4AllTD.Coord.perSens.x_null{i_win} = x_null';
    Kprime4AllTD.Coord.perSens.y_null{i_win} = y_null';
    Kprime4AllTD.Coord.perSens.z_null{i_win} = z_null';
end

%% 3) Compute Vector Norm and Vector Angle for each sensor and for predetermined sensor-clusters
vector_InformationLine = param.KprimeComparison.coord_EndPoint_InformationLine;
vector_DurationLine = param.KprimeComparison.coord_EndPoint_DurationLine;

for i_win = 1:size(Kprime4AllTD.Coord.perSens.x_exp,2) %3 Time windows
    %Experimental data
    %Accumulate coordinates in vector
    for i_sens = 1:size(Kprime4AllTD.Coord.perSens.x_exp{i_win},2)
        Kprime4AllTD.Vector.perSens.exp{i_win}(i_sens,:) = ...
            [Kprime4AllTD.Coord.perSens.x_exp{i_win}(i_sens) ...
            Kprime4AllTD.Coord.perSens.y_exp{i_win}(i_sens) ...
            Kprime4AllTD.Coord.perSens.z_exp{i_win}(i_sens)];
        
        %Compute vector norm (Euclidean length) of the vector
        %connecting origin (0,0,0) and data point (k' for TD: 150,300,600)
        Kprime4AllTD.Norm.perSens.exp{i_win}(i_sens,1) = ...
            norm(Kprime4AllTD.Vector.perSens.exp{i_win}(i_sens,:));
        
        %Compute vector angle to Information line
        Kprime4AllTD.AngleInformationLine.perSens.exp{i_win}(i_sens,1) = ...
            atan2(norm(cross(Kprime4AllTD.Vector.perSens.exp{i_win}(i_sens,:),vector_InformationLine)),...
            dot(Kprime4AllTD.Vector.perSens.exp{i_win}(i_sens,:),vector_InformationLine));
        
        %Compute vector angle to Duration line
        Kprime4AllTD.AngleDurationLine.perSens.exp{i_win}(i_sens,1) = ...
            atan2(norm(cross(Kprime4AllTD.Vector.perSens.exp{i_win}(i_sens,:),vector_DurationLine)),...
            dot(Kprime4AllTD.Vector.perSens.exp{i_win}(i_sens,:),vector_DurationLine));
        
        %Shuffled data
        for i_rep = 1:size(Kprime4AllTD.Coord.perSens.x_null{1},1)
            Kprime4AllTD.Vector.perSens.shuff{i_win}(i_sens,:,i_rep) = ...
                [Kprime4AllTD.Coord.perSens.x_null{i_win}(i_rep,i_sens) ...
                Kprime4AllTD.Coord.perSens.y_null{i_win}(i_rep,i_sens) ...
                Kprime4AllTD.Coord.perSens.z_null{i_win}(i_rep,i_sens)];
            
            %Compute vector norm (Euclidean length) of the vector
            %connecting origin (0,0,0) and data point (k' for TD: 150,300,600)
            Kprime4AllTD.Norm.perSens.shuff{i_win}(i_sens,i_rep) = ...
                norm(Kprime4AllTD.Vector.perSens.shuff{i_win}(i_sens,:,i_rep));
            
            %Compute vector angle to Information line
            Kprime4AllTD.AngleInformationLine.perSens.shuff{i_win}(i_sens,i_rep) = ...
                atan2(norm(cross(Kprime4AllTD.Vector.perSens.shuff{i_win}(i_sens,:,i_rep),vector_InformationLine)),...
                dot(Kprime4AllTD.Vector.perSens.shuff{i_win}(i_sens,:,i_rep),vector_InformationLine));
            
            %Compute vector angle to Duration line
            Kprime4AllTD.AngleDurationLine.perSens.shuff{i_win}(i_sens,i_rep) = ...
                atan2(norm(cross(Kprime4AllTD.Vector.perSens.shuff{i_win}(i_sens,:,i_rep),vector_DurationLine)),...
                dot(Kprime4AllTD.Vector.perSens.shuff{i_win}(i_sens,:,i_rep),vector_DurationLine));
        end
    end
end

%Average vector norm and angle across sensors in specific cluster
%Load in file containing all sensor-clusters
load([paths_NASTD_MEG.Current_analysis 'clusterSOI.mat'])

for i_win = 1:size(Kprime4AllTD.Coord.perSens.x_exp,2) %3 Time windows
    for i_clustertype = 1:length(param.KprimeComparison.ClusterSOI)
        clusterlabel = ['per' param.KprimeComparison.ClusterSOI{i_clustertype}];
        for i_clusterwindow = 1:length(ClusterSOI.(param.KprimeComparison.ClusterSOI{i_clustertype}))
            for i_clusternum = 1:length(ClusterSOI.(param.KprimeComparison.ClusterSOI{i_clustertype}){i_clusterwindow})
                if ~isempty(ClusterSOI.(param.KprimeComparison.ClusterSOI{i_clustertype})...
                        {i_clusterwindow}{i_clusternum}')
                    %Experimental data
                    %Vector norm
                    Kprime4AllTD.Norm.(clusterlabel).exp{i_clusterwindow}{i_clusternum} = ...
                        mean(Kprime4AllTD.Norm.perSens.exp{i_clusterwindow}...
                        (ClusterSOI.(param.KprimeComparison.ClusterSOI{i_clustertype})...
                        {i_clusterwindow}{i_clusternum}'));
                    %vector angle to Information line
                    Kprime4AllTD.AngleInformationLine.(clusterlabel).exp{i_clusterwindow}{i_clusternum} = ...
                        mean(Kprime4AllTD.AngleInformationLine.perSens.exp{i_clusterwindow}...
                        (ClusterSOI.(param.KprimeComparison.ClusterSOI{i_clustertype})...
                        {i_clusterwindow}{i_clusternum}'));
                    %vector angle to Duration line
                    Kprime4AllTD.AngleDurationLine.(clusterlabel).exp{i_clusterwindow}{i_clusternum} = ...
                        mean(Kprime4AllTD.AngleDurationLine.perSens.exp{i_clusterwindow}...
                        (ClusterSOI.(param.KprimeComparison.ClusterSOI{i_clustertype})...
                        {i_clusterwindow}{i_clusternum}'));
                    
                    %Shuffled data
                    %Vector norm
                    Kprime4AllTD.Norm.(clusterlabel).shuff{i_clusterwindow}{i_clusternum} = ...
                        mean(Kprime4AllTD.Norm.perSens.shuff{i_clusterwindow}...
                        (ClusterSOI.(param.KprimeComparison.ClusterSOI{i_clustertype})...
                        {i_clusterwindow}{i_clusternum}',:));
                    %vector angle to Information line
                    Kprime4AllTD.AngleInformationLine.(clusterlabel).shuff{i_clusterwindow}{i_clusternum} = ...
                        mean(Kprime4AllTD.AngleInformationLine.perSens.shuff{i_clusterwindow}...
                        (ClusterSOI.(param.KprimeComparison.ClusterSOI{i_clustertype})...
                        {i_clusterwindow}{i_clusternum}',:));
                    %vector angle to Duration line
                    Kprime4AllTD.AngleDurationLine.(clusterlabel).shuff{i_clusterwindow}{i_clusternum} = ...
                        mean(Kprime4AllTD.AngleDurationLine.perSens.shuff{i_clusterwindow}...
                        (ClusterSOI.(param.KprimeComparison.ClusterSOI{i_clustertype})...
                        {i_clusterwindow}{i_clusternum}',:));
                end
            end
        end
    end
end

%Compute nonparametrical p-values for Norm and Angle for experimental data
for i_win = 1:size(Kprime4AllTD.Coord.perSens.x_exp,2) %3 Time windows
    for i_sens = 1:size(Kprime4AllTD.Coord.perSens.x_exp{i_win},2)
        %per sensor
        %Vector norm (p-val determined by shuff > exp)
        Kprime4AllTD.Norm.perSens.pval_expvsshuff{i_win}(i_sens,1) = ...
            sum(Kprime4AllTD.Norm.perSens.shuff{i_win}(i_sens,:) >= ...
            (Kprime4AllTD.Norm.perSens.exp{i_win}(i_sens))) ...
            / size(Kprime4AllTD.Coord.perSens.x_null{1},1);
        %Vector angle to Information line (p-val determined by shuff < exp)
        Kprime4AllTD.AngleInformationLine.perSens.pval_expvsshuff{i_win}(i_sens,1) = ...
            sum(Kprime4AllTD.AngleInformationLine.perSens.shuff{i_win}(i_sens,:) <= ...
            (Kprime4AllTD.AngleInformationLine.perSens.exp{i_win}(i_sens))) ...
            / size(Kprime4AllTD.Coord.perSens.x_null{1},1);
        %Vector angle to Duration line (p-val determined by shuff < exp)
        Kprime4AllTD.AngleDurationLine.perSens.pval_expvsshuff{i_win}(i_sens,1) = ...
            sum(Kprime4AllTD.AngleDurationLine.perSens.shuff{i_win}(i_sens,:) <= ...
            (Kprime4AllTD.AngleDurationLine.perSens.exp{i_win}(i_sens))) ...
            / size(Kprime4AllTD.Coord.perSens.x_null{1},1);
    end
    %per sensor-cluster
    for i_clustertype = 1:length(param.KprimeComparison.ClusterSOI)
        clusterlabel = ['per' param.KprimeComparison.ClusterSOI{i_clustertype}];
        for i_clusterwindow = 1:length(ClusterSOI.(param.KprimeComparison.ClusterSOI{i_clustertype}))
            for i_clusternum = 1:length(ClusterSOI.(param.KprimeComparison.ClusterSOI{i_clustertype}){i_clusterwindow})
                if ~isempty(ClusterSOI.(param.KprimeComparison.ClusterSOI{i_clustertype})...
                        {i_clusterwindow}{i_clusternum}')
                    %Experimental data
                    %Vector norm (p-val determined by shuff > exp)
                    Kprime4AllTD.Norm.(clusterlabel).pval_expvsshuff{i_clusterwindow}{i_clusternum} = ...
                        sum(Kprime4AllTD.Norm.(clusterlabel).shuff{i_clusterwindow}{i_clusternum} >= ...
                        (Kprime4AllTD.Norm.(clusterlabel).exp{i_clusterwindow}{i_clusternum})) ...
                        / size(Kprime4AllTD.Coord.perSens.x_null{1},1);                  
                    %vector angle to Information line (p-val determined by shuff < exp)
                    Kprime4AllTD.AngleInformationLine.(clusterlabel).pval_expvsshuff{i_clusterwindow}{i_clusternum} = ...
                        sum(Kprime4AllTD.AngleInformationLine.(clusterlabel).shuff{i_clusterwindow}{i_clusternum} <= ...
                        (Kprime4AllTD.AngleInformationLine.(clusterlabel).exp{i_clusterwindow}{i_clusternum})) ...
                        / size(Kprime4AllTD.Coord.perSens.x_null{1},1);
                    %vector angle to Duration line (p-val determined by shuff < exp)
                    Kprime4AllTD.AngleDurationLine.(clusterlabel).pval_expvsshuff{i_clusterwindow}{i_clusternum} = ...
                        sum(Kprime4AllTD.AngleDurationLine.(clusterlabel).shuff{i_clusterwindow}{i_clusternum} <= ...
                        (Kprime4AllTD.AngleDurationLine.(clusterlabel).exp{i_clusterwindow}{i_clusternum})) ...
                        / size(Kprime4AllTD.Coord.perSens.x_null{1},1);
                end
            end
        end
    end
end

%% 4) Save output data
savefile = [path_outputdata 'Group_KprimeNormAngle.mat'];
save(savefile, ...
    'Kprime4AllTD', ... 
    '-v7.3');
end
