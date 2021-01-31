function NASTD_MEG_SHI_PlotKprimeNormAngle_Group...
    (param, ...
    paths_NASTD_MEG)
%Aim: Plot ouput of statistical exp vs. shuff comparison for exp vs. shuff vector norm and angle
%Plot histograms for exp vs shuff Vector Norm and Angle for all predetermined sensor cluster
%Plot 3D plot for experimental vs shuffled data
%% 0) Specify vars, paths, and setup fieldtrip
addpath(genpath(paths_NASTD_MEG.ScriptsDir));

%Data output
path_inputdata = [paths_NASTD_MEG.Current_outputdata 'ExpK/TDcomparison/'];

%Load in predetermined cluster-sensor indices
load([paths_NASTD_MEG.Current_analysis 'clusterSOI.mat'])

%% 1) Load in data
%Data includes Kprime coordinates, vector norm and angle data, and exp vs.
%shuff p-values per sensor and predetermiend sensor-clusters
load([path_inputdata 'Group_KprimeNormAngle.mat']);

%% 2) Determine analysis parameters
%Define time windows used for analysis
win_size    = 30; 
win_overlap = 0; 
samplingFreq = 600;
nSamplesPerTone = 0.15 * samplingFreq; %shortest tone dur * samplingFreq

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

%% 3) Plot histogram showing exp vs. shuff data per cluster
for i_clustertype = 1:length(param.KprimeComparison.ClusterSOI)
    curr_clustertyp = param.KprimeComparison.ClusterSOI{i_clustertype};
    label_curr_clustertyp = ['per' curr_clustertyp];
    disp(['Plotting histograms for Kprime data computed for cluster: ' curr_clustertyp])
    
    %Determine number of clusters per time window
    for i_win = 1:length(ClusterSOI.(curr_clustertyp))
        maxNumCluster_perTW(i_win) = length(ClusterSOI.(curr_clustertyp){i_win});
    end
    maxNumCluster = max(maxNumCluster_perTW);
    
    for i_win = 1:length(ClusterSOI.(curr_clustertyp))
        
        %Determine TW parameters for text plotting
        w1 = num2str(1000 * windows(i_win, 1) / samplingFreq , 3 );
        w2 = num2str(1000 * windows(i_win, 2) / samplingFreq , 3 );
        win_text = [w1 '-' w2 'ms'];
        
        %Plot histograms for vector norm, vector angle to information line, vector angle to duration line
        if ~isempty(ClusterSOI.(curr_clustertyp){i_win}{1})
            
            h = figure;
            set(gcf,'units','normalized','outerposition',[0 0 1 1])
            subplot_position = 0;
            
            for i_cluster = 1:length(ClusterSOI.(curr_clustertyp){i_win})
                indices_sensorcluster = ClusterSOI.(param.KprimeComparison.ClusterSOI{i_clustertype})...
                    {i_win}{i_cluster}';
                
                if ~isempty(indices_sensorcluster)
                    %Vector Norm
                    subplot_position = subplot_position + 1;
                    subplot(maxNumCluster,3,subplot_position);
                    %plot histogram for shuffled values
                    hist(mean(Kprime4AllTD.Norm.perSens.shuff{i_win}(indices_sensorcluster,:),1),20);
                    h = findobj(gca,'Type','patch');
                    h.FaceColor = [0.5 0.5 0.5]; h.EdgeColor = 'k'; hold on;
                    %plot experimental value
                    plot([mean(Kprime4AllTD.Norm.perSens.exp{i_win}(indices_sensorcluster)) ...
                        mean(Kprime4AllTD.Norm.perSens.exp{i_win}(indices_sensorcluster))],[0 175],'g','LineWidth',5);
                    title(['Vector Norm - Cluster ' num2str(i_cluster)]);
                    legend(['p = ' num2str(Kprime4AllTD.Norm.(label_curr_clustertyp).pval_expvsshuff{i_win}{i_cluster})])
                    xlabel('Vector Norm');
                    xlim([7 12])
                    set(gca,'XTick', [7 8 9 10 11 12])
                    ylabel(['Freq (' num2str(length(Kprime4AllTD.Norm.perSens.shuff{i_win})) ' permutations)'])
                    ylim([0 175])
                    set(gca,'YTick', [0 50 100 150])
                    
                    %Vector Angle to Information Line
                    subplot_position = subplot_position + 1;
                    subplot(maxNumCluster,3,subplot_position);
                    %plot histogram for shuffled values
                    hist(mean(Kprime4AllTD.AngleInformationLine.perSens.shuff{i_win}(indices_sensorcluster,:),1),20);
                    h = findobj(gca,'Type','patch');
                    h.FaceColor = [0.5 0.5 0.5]; h.EdgeColor = 'k'; hold on;
                    %plot experimental value
                    plot([mean(Kprime4AllTD.AngleInformationLine.perSens.exp{i_win}(indices_sensorcluster)) ...
                        mean(Kprime4AllTD.AngleInformationLine.perSens.exp{i_win}(indices_sensorcluster))],[0 175],'g','LineWidth',5);
                    title(['Vector Angle to Information Line - Cluster ' num2str(i_cluster)]);
                    legend(['p = ' num2str(Kprime4AllTD.AngleInformationLine.(label_curr_clustertyp).pval_expvsshuff{i_win}{i_cluster})])
                    xlabel('Vector Angle [radians]');
                    xlim([0 0.2])
                    set(gca,'XTick', [0 0.05 0.1 0.15 0.2])
                    ylabel(['Freq (' num2str(length(Kprime4AllTD.Norm.perSens.shuff{i_win})) ' permutations)'])
                    ylim([0 175])
                    set(gca,'YTick', [0 100 150])
                    
                    subplot_position = subplot_position + 1;
                    subplot(maxNumCluster,3,subplot_position);
                    %plot histogram for shuffled values
                    hist(mean(Kprime4AllTD.AngleDurationLine.perSens.shuff{i_win}(indices_sensorcluster,:),1),20);
                    h = findobj(gca,'Type','patch');
                    h.FaceColor = [0.5 0.5 0.5]; h.EdgeColor = 'k'; hold on;
                    %plot experimental value
                    plot([mean(Kprime4AllTD.AngleDurationLine.perSens.exp{i_win}(indices_sensorcluster)) ...
                        mean(Kprime4AllTD.AngleDurationLine.perSens.exp{i_win}(indices_sensorcluster))],[0 175],'g','LineWidth',5);
                    title(['Vector Angle to Duration Line - Cluster ' num2str(i_cluster)]);
                    legend(['p = ' num2str(Kprime4AllTD.AngleDurationLine.(label_curr_clustertyp).pval_expvsshuff{i_win}{i_cluster})])
                    xlabel('Vector Angle [radians]');
                    xlim([0.3 0.6])
                    set(gca,'XTick', [0.3 0.4 0.5 0.6])
                    ylabel(['Freq (' num2str(length(Kprime4AllTD.Norm.perSens.shuff{i_win})) ' permutations)'])
                    ylim([0 175])
                    set(gca,'YTick', [0 100 150])
                end
            end
            %Add overall title
            t = suptitle(['Vector Norm + Angle for Kprime values (exp vs. shuff data) computed for ' curr_clustertyp '; TW: ' win_text]);

            t.Interpreter = 'none';
            
            if param.plot.save   == 1
                path_figs = [paths_NASTD_MEG.Current_outputfig 'ExpKprime/TDcomparison/Histogram/'];
                mkdir(path_figs)
                
                filename     = ['Group_HistogramKprime_' curr_clustertyp '_' win_text '.png'];
                figfile      = [path_figs filename];
                saveas(gcf, [figfile], 'png'); %save png version
                close;
            end
        end
    end
end

%% 4) Plot 3D representation of exp vs. shuff data per cluster
for i_clustertype = 1:length(param.KprimeComparison.ClusterSOI)
    curr_clustertyp = param.KprimeComparison.ClusterSOI{i_clustertype};
    label_curr_clustertyp = ['per' curr_clustertyp];
    disp(['Plotting 3D representations for Kprime data computed for cluster: ' curr_clustertyp])
    
    vector_InformationLine = param.KprimeComparison.coord_EndPoint_InformationLine;
    vector_DurationLine = param.KprimeComparison.coord_EndPoint_DurationLine;
    
    for i_win = 1:length(ClusterSOI.(curr_clustertyp))
        
        %Determine TW parameters for text plotting
        w1 = num2str(1000 * windows(i_win, 1) / samplingFreq , 3 );
        w2 = num2str(1000 * windows(i_win, 2) / samplingFreq , 3 );
        win_text = [w1 '-' w2 'ms'];
        
        for i_cluster = 1:length(ClusterSOI.(param.KprimeComparison.ClusterSOI{i_clustertype}){i_win})
            
            indices_sensorcluster = ClusterSOI.(param.KprimeComparison.ClusterSOI{i_clustertype})...
                {i_win}{i_cluster}';
            
            if ~isempty(indices_sensorcluster)
                %Plot 3D representation of exp vs. shuff data per cluster
                h = figure;
                set(gcf,'units','normalized','outerposition',[0 0 1 1])
                set(gcf,'Renderer','painters');
                
                %Plot hypothesis-derived orientation lines
                plot3([0:vector_InformationLine(1)],[0:vector_InformationLine(2)],[0:vector_InformationLine(3)], ...
                    '-', 'color', [0.4 0.8 1.0],'LineWidth', 5) %Information bottleneck: "Information Line" (Hypothesis 2)
                hold on
                plot3([0:1:vector_DurationLine(1)],[0:0.5:vector_DurationLine(2)],[0:0.25:vector_DurationLine(3)], ...
                    '-', 'color', [1.0 0.4 0.4],'LineWidth', 5) %Temporal bottleneck: "Duration Line" (Hypothesis 1)
                grid on
                xlabel('Kprime (150 ms TD)')
                ylabel('Kprime (300 ms TD)')
                zlabel('Kprime (600 ms TD)')
                set(gca,'XTick', [0 2 4 6 8])
                set(gca,'YTick', [0 2 4 6 8])
                set(gca,'ZTick', [0 2 4 6 8])
                set(gca,'Xlim', [0 8])
                set(gca,'Ylim', [0 8])
                set(gca,'Zlim', [0 8])
                axis square
                set(gca,'linewidth',3)
                view(-15,15)
                
                %Plot center of mass (across sensors) for shuffled Kprime values for each repetition
                scatter3(...
                    mean(Kprime4AllTD.Coord.perSens.x_null{i_win}(:,indices_sensorcluster),2), ...
                    mean(Kprime4AllTD.Coord.perSens.y_null{i_win}(:,indices_sensorcluster),2), ...
                    mean(Kprime4AllTD.Coord.perSens.z_null{i_win}(:,indices_sensorcluster),2), ...
                    45,  [0.2 0.2 0.2], 'filled', ...
                    'd', 'MarkerFaceAlpha', 0.5);
                %Plot border around shuffled Kprime value distribution
                Xcoord_allShuff = [];
                Ycoord_allShuff = [];
                Zcoord_allShuff = [];
                for i_repetition = 1:length(Kprime4AllTD.Coord.perSens.x_null{i_win})
                    Xcoord_allShuff = [Xcoord_allShuff; Kprime4AllTD.Coord.perSens.x_null{i_win}(i_repetition, indices_sensorcluster)'];
                    Ycoord_allShuff = [Ycoord_allShuff; Kprime4AllTD.Coord.perSens.y_null{i_win}(i_repetition, indices_sensorcluster)'];
                    Zcoord_allShuff = [Zcoord_allShuff; Kprime4AllTD.Coord.perSens.z_null{i_win}(i_repetition, indices_sensorcluster)'];
                end
                border_pointcloudSHUFFclusterSOI = ...
                    boundary([Xcoord_allShuff, Ycoord_allShuff, Zcoord_allShuff]);
                hold on;
                trisurf(border_pointcloudSHUFFclusterSOI, ...
                    Xcoord_allShuff, ...
                    Ycoord_allShuff, ...
                    Zcoord_allShuff, ...
                    'FaceColor', [0.2 0.2 0.2], 'FaceAlpha', 0.1, 'EdgeAlpha', 0.1);
                
                %Plot experimental Kprime values for respective selected sensors
                scatter3(...
                    Kprime4AllTD.Coord.perSens.x_exp{i_win}(indices_sensorcluster), ...
                    Kprime4AllTD.Coord.perSens.y_exp{i_win}(indices_sensorcluster), ...
                    Kprime4AllTD.Coord.perSens.z_exp{i_win}(indices_sensorcluster), ...
                    45,  [0.109 0.596, 0.086], 'filled', ...
                    'MarkerFaceAlpha', 0.5);
                %Plot center of mass across sensors for experimental Kprime distribution
                scatter3(...
                    mean(Kprime4AllTD.Coord.perSens.x_exp{i_win}(indices_sensorcluster)), ...
                    mean(Kprime4AllTD.Coord.perSens.y_exp{i_win}(indices_sensorcluster)), ...
                    mean(Kprime4AllTD.Coord.perSens.z_exp{i_win}(indices_sensorcluster)), ...
                    250,  [0.109 0.596, 0.086], 'filled', ...
                    'd', 'MarkerFaceAlpha', 1);
                %Plot border around around experimental Kprime value
                border_pointcloudEXPclusterSOI = ...
                    boundary([Kprime4AllTD.Coord.perSens.x_exp{i_win}(indices_sensorcluster); ...
                    Kprime4AllTD.Coord.perSens.y_exp{i_win}(indices_sensorcluster); ...
                    Kprime4AllTD.Coord.perSens.z_exp{i_win}(indices_sensorcluster)]');
                hold on;
                trisurf(border_pointcloudEXPclusterSOI, ...
                    Kprime4AllTD.Coord.perSens.x_exp{i_win}(indices_sensorcluster), ...
                    Kprime4AllTD.Coord.perSens.y_exp{i_win}(indices_sensorcluster), ...
                    Kprime4AllTD.Coord.perSens.z_exp{i_win}(indices_sensorcluster), ...
                    'FaceColor', [0.109 0.596, 0.086], 'FaceAlpha', 0.1, 'EdgeAlpha', 0.1)
                
                title({['Kprime values (exp vs. shuff data) computed for ' curr_clustertyp '(Cluster: ' num2str(i_cluster) '); TW: ' win_text],...
                    ['Norm: p = ' num2str(Kprime4AllTD.Norm.(label_curr_clustertyp).pval_expvsshuff{i_win}{i_cluster}), ...
                    '; Angle (Information Line): p = ' num2str(Kprime4AllTD.AngleInformationLine.(label_curr_clustertyp).pval_expvsshuff{i_win}{i_cluster}), ...
                    '; Angle (Duration Line): p = ' num2str(Kprime4AllTD.AngleDurationLine.(label_curr_clustertyp).pval_expvsshuff{i_win}{i_cluster})]}, ...
                    'Interpreter','none')
                
                if param.plot.save   == 1
                    path_figs = [paths_NASTD_MEG.Current_outputfig 'ExpKprime/TDcomparison/3Drepresentation/'];
                    mkdir(path_figs)
                    
                    filename     = ['Group_3DKprime_' curr_clustertyp '_Cluster' num2str(i_cluster) '_' win_text '.png'];
                    figfile      = [path_figs filename];
                    saveas(gcf, [figfile], 'png'); %save png version
                    delete(h);
                end
            end
        end
    end
end

end