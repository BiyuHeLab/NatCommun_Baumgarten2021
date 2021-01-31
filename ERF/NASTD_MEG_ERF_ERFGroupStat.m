function GroupStat_ERF = NASTD_MEG_ERF_ERFGroupStat...
    (i_tonedur, param, paths_NASTD_MEG)

%Aim: Compute group-level statistical comparison between low vs. high p*34 trials

%Load input data and aggregate across subjects
path_input = [paths_NASTD_MEG.Current_outputdata 'Data_PredProcCluster/Ssub/' ...
    param.tonedur_text{i_tonedur} 'sTD/'];

ERF_avgTrialsperCond_Allsub = struct;

for i_sub = 1:length(param.subs)

 load([path_input  param.subs{i_sub} '_PredProcCluster_' ...
        param.tonedur_text{i_tonedur} 'sTD.mat']);
    
    conditions = fieldnames(ERF_avgTrialsperCond);
    for i_cond = 1:length(conditions)
        ERF_avgTrialsperCond_Allsub.(conditions{i_cond}){i_sub} = ...
            ERF_avgTrialsperCond.(conditions{i_cond});
    end
end

%Determine samples corresponding to the respective tone starts
nSamplesPerTone = str2double(param.tonedur_text{i_tonedur}) * param.samplefreq;
index_T1start = find(round(ERF_avgTrialsperCond_Allsub.predp34_low{1}.time,4) == 0);
Samples_AllToneStart = [index_T1start : nSamplesPerTone : (index_T1start + (34*nSamplesPerTone))];

conditions = fieldnames(ERF_avgTrialsperCond_Allsub);
for i_cond = 1:length(conditions)   
  
    %Compute grand average for plotting
    cfg                 = [];
    cfg.parameter       = 'avg';
    cfg.keepindividual  = 'no ';
    GroupStat_ERF.avg.(conditions{i_cond}) = ...
        ft_timelockgrandaverage(cfg, ERF_avgTrialsperCond_Allsub.(conditions{i_cond}){:});    
    
end

%Determine analysis window (start first tone to start response window)
TOI_statcomp = ...
    [GroupStat_ERF.avg.predp34_low.time(Samples_AllToneStart(1))...%start first tone
    GroupStat_ERF.avg.predp34_low.time(Samples_AllToneStart(35)) + 0.4]; %end last tone + 400ms / start response window

%Compute Group-level statistics (nonparametric, cluster corrected)
cfg                     = [];
cfg.channel             = 'all'; %cluster-SOI have been preselected
cfg.latency             = TOI_statcomp; %start first tone to start response window
cfg.avgovertime         = 'no';
cfg.parameter           = 'avg';

cfg.method              = 'montecarlo';
cfg.alpha               = param.stat.pval_clusterstat;
cfg.tail                = 0;        %0 = two-sided test
cfg.correcttail         = 'prob';   %Correct prob subfield to reflect two-sided test

cfg.correctm            = 'cluster';
cfg.clusteralpha        = 0.05;
cfg.clusterstatistic    = 'maxsum';
cfg.clustertail         = 0;

cfg.numrandomization    = param.stat.NumReps;

Nsub                    = length(param.subs);
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];

cfg.statistic           = 'ft_statfun_depsamplesT'; %within-subjects comparison

GroupStat_ERF.stat_LowvsHighPredp34 = ft_timelockstatistics(cfg,...
    ERF_avgTrialsperCond_Allsub.predp34_low{:},...
    ERF_avgTrialsperCond_Allsub.predp34_high{:});

end