function NASTD_MEG_Corr_PlotCorrSOI ...
    (ExpKprime_AllSub_AvgToneDur, Fstatinteraction_FTPLrating, ...
    ClusterSOI, ...
    rho_CorrexpKFstat_perSOI, pval_CorrexpKFstat_perSOI, pvalFDRcorrected_CorrexpKFstat_perSOI, ...
    param_NASTD_MEG, ...
    paths_NASTD_MEG)

%Aim: Plot correlation between Kprime values (averaged across clusterSOI)
%and FTPLrating (Fvalues) across subjects in 2D subplots (1 plot per TD,
%with TW as rows and Clusters as columns

%% 1) Determine parameters
%Determine Window length
win_size    = 30;
win_overlap = 0;
samplingFreq = 600;
toneDur_inSecs = 0.15;
nSamplesPerTone = toneDur_inSecs * samplingFreq;
windows = [1 win_size];

while windows(end,end) < nSamplesPerTone
    windows = [windows; windows(end,:) + (win_size - win_overlap)];
end
if windows(end,end) > nSamplesPerTone
    windows(end,:) = [];
end
windows_inms = (windows / samplingFreq) * 1000;

for i_clustertype = 1:length(param_NASTD_MEG.KprimeComparison.ClusterSOI)
    
    clusterlabel = [param_NASTD_MEG.KprimeComparison.ClusterSOI{i_clustertype}];
    
    %Find max. number of cluster per TW (to determine number rows in subplot)
    for i_win = 1:length(ClusterSOI.(clusterlabel)) %Loop across timewin (lowest common window length)
        maxNumCluster_perTW(i_win) = length(ClusterSOI.(clusterlabel){i_win});
    end
    maxNumCluster = max(maxNumCluster_perTW);
    
    %Find limits for axes
    Xlim = 0;
    for i_win = 1:length(rho_CorrexpKFstat_perSOI)
        for i_cluster = 1:length(ClusterSOI.(clusterlabel){i_win})
            if ceil(max(mean(ExpKprime_AllSub_AvgToneDur{i_win}(:,ClusterSOI.(clusterlabel){i_win}{i_cluster}),2))) > Xlim
                Xlim = ceil(max(mean(ExpKprime_AllSub_AvgToneDur{i_win}(:,ClusterSOI.(clusterlabel){i_win}{i_cluster}),2)));
            end
        end
    end
    Ylim = ceil(max(Fstatinteraction_FTPLrating));
    
    %% 2) Set up figure struct
    h = figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    subplot_position = 0;
    
    for i_win = 1:length(rho_CorrexpKFstat_perSOI)
        
        %Prepare struct for plot and plot
        w1 = num2str( 1000 * windows(i_win, 1) / samplingFreq , 3 );
        w2 = num2str( 1000 * windows(i_win, 2) / samplingFreq , 3 );
        win_title = ['TW = [' w1 '-' w2 'ms]'];
        
        for i_cluster = 1:length(ClusterSOI.(clusterlabel){i_win})
            
            %Read out data points
            plotK = mean(ExpKprime_AllSub_AvgToneDur{i_win}(:,ClusterSOI.(clusterlabel){i_win}{i_cluster}),2);
            plotF = Fstatinteraction_FTPLrating;
            
            %compute polynomial fit
            fitPts = polyfit(plotK,plotF,1);
            minx = 0;
            maxx = max(plotK)+1;
            miny = minx*fitPts(1) + fitPts(2);
            maxy = maxx*fitPts(1) + fitPts(2);
            
            %set subplot position
            subplot_position = subplot_position + 1;
            subplot(length(rho_CorrexpKFstat_perSOI), maxNumCluster, subplot_position);
            
            %plot data points
            scatter(plotK, plotF, 150, 'filled', 'k', 'MarkerFaceAlpha', 0.6);
            hold on
            plot([minx maxx], [miny, maxy])
%             xlim([0 Xlim])
%             ylim([0, Ylim])
            xlabel('K')
            ylabel('F')
            
            %Add subtitle
            t = title({[win_title '; Cluster#: ' num2str(i_cluster)],...
                ['r = ' num2str(round(rho_CorrexpKFstat_perSOI{i_win,i_cluster},2)) ...
                '; p = ' num2str(round(pval_CorrexpKFstat_perSOI{i_win,i_cluster},3)) ...
                '; p(FDR) =' num2str(round(pvalFDRcorrected_CorrexpKFstat_perSOI{i_win,i_cluster},3))]});
            set(t,'FontSize',12)
            
            %Add figure title
            s = suptitle(['Correlation Kprime (' clusterlabel ...
                ') & Fstatistic for behavioral p34 x p*34 interaction effect across subjects']);
            set(s,'FontSize',16)
            set(s,'Interpreter', 'none')
            
            %Resize
            pos = get(gca, 'Position');
            pos(3) = pos(3)*0.75;
            pos(4) = pos(4)*0.75;
            set(gca, 'Position', pos);
        end
    end
    
    %% 3) Save figure
    %Save figure
    if param_NASTD_MEG.plot.save   == 1
        path_figs = [paths_NASTD_MEG.Current_outputfig '2Dplot_clusterSOI/'];
        mkdir(path_figs)
        
        filename     = ['Group_2DCorrKprimeFstat_' clusterlabel '.png'];
        figfile      = [path_figs filename];
        saveas(gcf, [figfile], 'png'); %save png version
        delete(h);
    end
end

end