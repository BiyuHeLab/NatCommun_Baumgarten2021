function  NASTD_MEG_Behav_FTPLperSeqBeta_Group...
    (subs, tonedur_text, ...
    paths_NASTD_MEG)

%Plot and Compute group level behavioral Final Tone Pitch Likelihood (FTPL) 
%rating depending on sequence beta levels across and for each tone duration

path_outputdata = ([paths_NASTD_MEG.Current_outputdata 'Statistics/']);
mkdir(path_outputdata);

seqBeta_text = {'0.5' '0.99' '1.5'};
p34_text = {'-3','-2','-1','1','2','3'};
predp34_text = {'low','med','high'};

%% 1) load behavioral data and prepare data for plotting
%Read out ftpl ratings per tone duration, p*34, p 34, and beta level
[subPred, subDurPred] = NASTD_MEG_Behav_FTPL_PrepData4ANOVA(subs, tonedur_text, paths_NASTD_MEG);

%Avg FTPL rating per seqBeta for each subject
for i_seqBeta = 1:length(seqBeta_text)
    for i_sub = 1:length(subs)
        
        %Across All Tone Durations
        AvgFTPLratingperSeqBeta{1}(i_sub,i_seqBeta) = ...
            nanmean(subPred{i_sub}.PredRating(strcmp(subPred{i_sub}.beta,seqBeta_text{i_seqBeta})));
        SDFTPLratingperSeqBeta{1}(i_sub,i_seqBeta) = ...
            nanstd(subPred{i_sub}.PredRating(strcmp(subPred{i_sub}.beta,seqBeta_text{i_seqBeta})));
        FTPLratingperSeqBeta{1}{i_sub}(:,i_seqBeta) = ...
            subPred{i_sub}.PredRating(strcmp(subPred{i_sub}.beta,seqBeta_text{i_seqBeta}));
               
        %Per Tone Duration
        for i_tonedur = 1:length(tonedur_text)-1
            AvgFTPLratingperSeqBeta{i_tonedur+1}(i_sub,i_seqBeta) = ...
                nanmean(subDurPred{i_sub,i_tonedur}.PredRating(strcmp(subDurPred{i_sub,i_tonedur}.beta,seqBeta_text{i_seqBeta})));
            SDFTPLratingperSeqBeta{i_tonedur+1}(i_sub,i_seqBeta) = ...
                nanstd(subDurPred{i_sub,i_tonedur}.PredRating(strcmp(subDurPred{i_sub,i_tonedur}.beta,seqBeta_text{i_seqBeta})));
            FTPLratingperSeqBeta{i_tonedur+1}{i_sub}(:,i_seqBeta) = ...
                subDurPred{i_sub,i_tonedur}.PredRating(strcmp(subDurPred{i_sub,i_tonedur}.beta,seqBeta_text{i_seqBeta}));
        end
    end
end

%Avg FTPL ratings across subjects
AvgFTPLratingperSeqBeta{1}(length(subs) + 1,:) = ...
    nanmean(AvgFTPLratingperSeqBeta{1});
SDFTPLratingperSeqBeta{1}(length(subs) + 1,:) = ...
    nanstd(AvgFTPLratingperSeqBeta{1}(1 : length(subs)));

for i_tonedur = 1:length(tonedur_text)-1
    AvgFTPLratingperSeqBeta{i_tonedur+1}(length(subs) + 1,:) = ...
        nanmean(AvgFTPLratingperSeqBeta{i_tonedur+1}(1 : length(subs));
    SDFTPLratingperSeqBeta{i_tonedur+1}(length(subs) + 1,:) = ...
        nanstd(AvgFTPLratingperSeqBeta{i_tonedur+1}(1 : length(subs)));
end

%% 2. PLot FTPL ratings per seq beta
%2.1 Plot FTPL ratings per seq beta across and for each tone duration (Fig. S7 aI)
h = figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
set(gcf,'Renderer','painters');

linecolor = {[0 0 0],[0.75, 0.75, 0.75],[0.25, 0.25, 0.25],[0.2, 0, 0]};

for i_tonedur = 1:length(tonedur_text)
    hold on;
    e = errorbar(1:3, [AvgFTPLratingperSeqBeta{i_tonedur}(length(subs) + 1,:)], ...
        [SDFTPLratingperSeqBeta{i_tonedur}(length(subs) + 1,:)/sqrt(length(subs))]);
    e.Color = linecolor{i_tonedur};
    e.LineWidth = 3;
    e.XData = e.XData+(i_tonedur*0.025);
end

ylim([2.75 3.75])
xlim([0.5 3.5])
set(gca,'XTick',[1 2 3])
set(gca,'XTickLabel', {'0.5','0.99','1.5'},'FontSize',12)
xlabel('Sequence Beta Level')
set(gca,'YTick',[2.75 3 3.25 3.5 3.75])
ylabel('FTPL rating','FontSize',12)
legend({'All TD', '150ms TD','300ms TD','600ms TD'})

suptitle(['Behavioral FTPL rating (Avg/SEM; n = ' num2str(length(subs)) ') for sequence Beta levels'])

%2.1 Plot FTPL ratings per p*34 for each seq beta, across and for each tone duration (Fig. S7 aII)
h = figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
set(gcf,'Renderer','painters');
CounterSubplot = 0;
DimSubplot = [3,1];

FTPLrating_perseqBeta_perSub_perp34predp34 = [];

for i_seqBeta = 1:length(seqBeta_text)
    
    FTPLrating = [];
    predp34 = [];
    p34 = [];
    subject = [];
    
    for i_sub = 1:length(subs)
        
        filter = strcmp(subPred{i_sub}.beta,seqBeta_text{i_seqBeta});
        
        FTPLrating = [FTPLrating; subPred{i_sub}.PredRating(filter)];
        predp34 = [predp34; subPred{i_sub}.predp34(filter)];
        p34 = [p34; subPred{i_sub}.p34(filter)];
        for i_trial = 1:length(subPred{i_sub}.PredRating(filter))
            subject = [subject; i_sub];
        end
        
        %Read out and average FTPLrating per seqBeta/predp34/p34 and sub
        for i_p34 = 1:length(unique(subPred{i_sub}.p34))
            for i_predp34 = 1:length(unique(subPred{i_sub}.predp34))
                
                filter = strcmp(subPred{i_sub}.beta,seqBeta_text{i_seqBeta}) & ...
                    strcmp(subPred{i_sub}.p34,p34_text{i_p34}) & ...
                    strcmp(subPred{i_sub}.predp34,predp34_text{i_predp34});
                
                FTPLrating_perseqBeta_perSub_perp34predp34{i_seqBeta}{i_predp34}(i_sub,i_p34) = ...
                    nanmean(subPred{i_sub}.PredRating(filter));
            end
        end
    end
    
    %Delete NaN entries (no response/missed trials)
    NAN_trials = find(isnan(FTPLrating));
    if ~isempty(NAN_trials)
        for i_NANtrials = NAN_trials
            FTPLrating(i_NANtrials) = [];
            predp34(i_NANtrials) = [];
            p34(i_NANtrials) = [];
            subject(i_NANtrials) = [];
        end
    end
    
    linecolor = {[0 0 1],[0 1 1],[1 1 0]};
    
    CounterSubplot = CounterSubplot + 1;
    sp = subplot(DimSubplot(1), DimSubplot(2),CounterSubplot);
    for i_predp34 = 1:length(predp34_text)
        hold on;
        e = errorbar(1:6, [mean(FTPLrating_perseqBeta_perSub_perp34predp34{i_seqBeta}{i_predp34})], ...
            [std(FTPLrating_perseqBeta_perSub_perp34predp34{i_seqBeta}{i_predp34})/sqrt(length(subs))]);
        e.Color = linecolor{i_predp34};
        e.LineWidth = 3;
    end
    ylim([1 5])
    xlim([0.5 6.5])
    set(gca,'XTick',[1 2 3 4 5 6])
    set(gca,'XTickLabel', {'220','277','349','554','698','880'},'FontSize',10)
    xlabel('p34 [Hz]')
    ylabel('FTPL rating across subjects','FontSize',8)
    title(['FTPLrating per p34 & p*34 for seqBeta: ' seqBeta_text{i_seqBeta}])
    if CounterSubplot == 1
        legend({'p*34: low','p*34: med','p*34: high'},'Location','best')
    end
    
end

suptitle(['Behavioral FTPLrating p34-p*34 interaction effect (Avg/SEM; n = ' ...
    num2str(length(subs)) ') per sequence Beta levels (Across TD)'])

%% 3) Group level FTPL rating comparisons via ANOVA: 
%Bring data in form for ANOVA
FTPLrating = [];
predp34 = [];
p34 = [];
seqBeta = [];
subject = [];
tonedur = [];

for i_sub = 1:length(subs)
    
    FTPLrating = [FTPLrating; subPred{i_sub}.PredRating];
    predp34 = [predp34; subPred{i_sub}.predp34_numerical];
    p34 = [p34; subPred{i_sub}.p34_numerical];
    seqBeta = [seqBeta; subPred{i_sub}.beta_numerical];
    tonedur = [tonedur; subPred{i_sub}.tonedur_numerical];

    for i_trial = 1:length(subPred{i_sub}.PredRating)
        subject = [subject; i_sub];
    end

end

%Delete NaN entries (no response/missed trials)
NAN_trials = find(isnan(FTPLrating));
if ~isempty(NAN_trials)
    for i_NANtrials = NAN_trials
        FTPLrating(i_NANtrials) = [];
        predp34(i_NANtrials) = [];
        p34(i_NANtrials) = [];
        seqBeta(i_NANtrials) = [];
        subject(i_NANtrials) = [];
        tonedur(i_NANtrials) = [];
    end
end

%Average FTPL ratings across trials for each subject
sumMat = [FTPLrating, p34, predp34, tonedur, seqBeta, subject];
SubAvgMat = nan(length(unique(p34))*length(unique(predp34))...
    *length(unique(tonedur))*length(unique(seqBeta))*length(subPred),6);
entry = 0;
for i_seqBeta = unique(seqBeta)'
    for i_Tonedur = unique(tonedur)'
        for i_predp34 = unique(predp34)'
            for i_p34 = unique(p34)'
                for i_sub = 1:length(subPred)
                    
                    entry = entry+1;
                    filter = sumMat(:,2) == i_p34 & ...
                        sumMat(:,3) == i_predp34 & ...
                        sumMat(:,4) == i_Tonedur & ...
                        sumMat(:,5) == i_seqBeta & ...
                        sumMat(:,6) == i_sub;
                    SubAvgMat(entry,1) = nanmean(sumMat(filter,1));
                    SubAvgMat(entry,2) = i_p34;
                    SubAvgMat(entry,3) = i_predp34;
                    SubAvgMat(entry,4) = i_Tonedur;
                    SubAvgMat(entry,5) = i_seqBeta;
                    SubAvgMat(entry,6) = i_sub;
                    
                end
            end
        end
    end
end

%Compute  3-way RM ANOVA (Factors: predp34, p34, ToneDur) for each beta level
for i_seqBeta = unique(seqBeta)'
    
    filter = SubAvgMat(:,5) == i_seqBeta;    
    [~,table_ANOVA3{i_seqBeta},~] = anovan(SubAvgMat(filter,1),...
        {SubAvgMat(filter,2),...
        SubAvgMat(filter,3),...
        SubAvgMat(filter,4),...
        SubAvgMat(filter,6)},...
        'model', 'full','random',4,'varnames',...
        {'p34','p*34','ToneDur','subjects'},'display','off');

end

end