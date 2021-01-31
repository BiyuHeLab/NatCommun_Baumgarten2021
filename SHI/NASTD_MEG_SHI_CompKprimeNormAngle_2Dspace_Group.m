function NASTD_MEG_SHI_CompKprimeNormAngle_2Dspace_Group...
    (param, ...
    paths_NASTD_MEG)
%Aim: Project 3D Kprime values on a 2D plane spanning between both
%hypothesis-derived orientation lines. Recompute p-values (exp vs
%shuff). Plot results in 2D coordinate system.

%% 0) Specify vars, paths, and setup fieldtrip
addpath(genpath(paths_NASTD_MEG.ScriptsDir));
%Data input
path_inputdata  = [paths_NASTD_MEG.Current_outputdata 'ExpK/TDcomparison/'];

%Data output
path_outputdata = [paths_NASTD_MEG.Current_outputdata 'ExpK/TDcomparison/2Dplane/'];
mkdir(path_outputdata)

%% 1) Load in data (Kprime coordinates )
load([path_inputdata 'Group_KprimeNormAngle.mat']);

%% 2) Construct 2D plane spanned by the orientation lines for number vs. duration
%Compute common plane
%Find equation (in standard form Ax+By+Cz = D) of plane passing through points
Point_null = param.KprimeComparison.coord_StartPoint;
Vector_InformationLine = param.KprimeComparison.coord_EndPoint_InformationLine;
Vector_DurationLine = param.KprimeComparison.coord_EndPoint_DurationLine;

% Compute normal vector to plane by computing cross-product
NormVectorPlane = cross(Vector_InformationLine,Vector_DurationLine); %Vector [A,B,C] defining plane

%Derive plane equation
syms x y z %construct symbolic vars
PointonPlane = [x, y, z];
plane_equation = dot(NormVectorPlane, PointonPlane-Point_null);
%compute dot product of vector between null and point on plane and normal vector
zplane = solve(plane_equation, z); %solve for z

% %Visual test if plane is aligned with orientation lines
% plot3([0:Vector_InformationLine(1)],[0:Vector_InformationLine(2)], [0:Vector_InformationLine(3)],...
%     '-', 'color', [0.4000    0.8000    1.0000], 'LineWidth', 5)
% hold on
% plot3([0:Vector_DurationLine(1)], [0:0.5:Vector_DurationLine(2)], [0:0.25:Vector_DurationLine(3)], '-',...
%     'color', [1.0000    0.4000    0.4000], 'LineWidth', 5)
% hold on
% fmesh(zplane)
% grid on
% xlabel('x - 150 ms TD')
% ylabel('y - 300 ms TD')
% zlabel('z - 600 ms TD')

%% 3) Project 3D Kprime points on 2D plane for all sensors
%Copy exp and shuff Kprime values to better format
for i_win = 1:size(Kprime4AllTD.Coord.perSens.x_exp,2) %3 Time windows
    for i_sensor = 1:size(Kprime4AllTD.Coord.perSens.x_exp{i_win},2)
        %Exp
        Kprime_exp_Original{i_win}(i_sensor,:) = ...
            [Kprime4AllTD.Coord.perSens.x_exp{i_win}(i_sensor) ...
            Kprime4AllTD.Coord.perSens.y_exp{i_win}(i_sensor) ...
            Kprime4AllTD.Coord.perSens.z_exp{i_win}(i_sensor)];
        
        %Shuff
        for i_rep = 1:length(Kprime4AllTD.Coord.perSens.x_null{i_win})
            Kprime_shuff_Original{i_win}{i_rep}(i_sensor,:) = ...
                [Kprime4AllTD.Coord.perSens.x_null{i_win}(i_rep,i_sensor) ...
                Kprime4AllTD.Coord.perSens.y_null{i_win}(i_rep,i_sensor) ...
                Kprime4AllTD.Coord.perSens.z_null{i_win}(i_rep,i_sensor)];
        end
    end
end

%Project exp and shuff Kprime values on plane
    %closestpoint: v = (d - sum(p.*n)) / sum(n.*n);
    %d = distance of plane from origin (0 since plane crosses point of origin)
    %p = point of interest
    %n = vector defining plane
for i_win = 1:length(Kprime4AllTD.Coord.perSens.x_exp)
    for i_sensor = 1:length(Kprime_exp_Original{i_win})
        
        %Exp
        closest_2dPoint = (0 - sum(Kprime_exp_Original{i_win}(i_sensor,:).*NormVectorPlane)) / sum(NormVectorPlane.*NormVectorPlane);
        KprimeEXP_OnPlane{i_win}(i_sensor,:) = Kprime_exp_Original{i_win}(i_sensor,:) + closest_2dPoint * NormVectorPlane;
        closest_2dPoint = [];
        
        Kprime4AllTD.Vector_2Dplane.perSens.exp{i_win}(i_sensor,1) = KprimeEXP_OnPlane{i_win}(i_sensor,1);
        Kprime4AllTD.Vector_2Dplane.perSens.exp{i_win}(i_sensor,2) = KprimeEXP_OnPlane{i_win}(i_sensor,2);
        Kprime4AllTD.Vector_2Dplane.perSens.exp{i_win}(i_sensor,3) = KprimeEXP_OnPlane{i_win}(i_sensor,3);
        
        %Shuff
        for i_rep = 1:length(Kprime4AllTD.Coord.perSens.x_null{i_win})
            closest_2dPoint = (0 - sum(Kprime_shuff_Original{i_win}{i_rep}(i_sensor,:).*NormVectorPlane)) / sum(NormVectorPlane.*NormVectorPlane);
            KprimeSHUFF_OnPlane{i_win}{i_rep}(i_sensor,:) = Kprime_shuff_Original{i_win}{i_rep}(i_sensor,:) + closest_2dPoint * NormVectorPlane;
            closest_2dPoint = [];
            
            Kprime4AllTD.Vector_2Dplane.perSens.shuff{i_win}(i_sensor,1,i_rep) = KprimeSHUFF_OnPlane{i_win}{i_rep}(i_sensor,1);
            Kprime4AllTD.Vector_2Dplane.perSens.shuff{i_win}(i_sensor,2,i_rep) = KprimeSHUFF_OnPlane{i_win}{i_rep}(i_sensor,2);
            Kprime4AllTD.Vector_2Dplane.perSens.shuff{i_win}(i_sensor,3,i_rep) = KprimeSHUFF_OnPlane{i_win}{i_rep}(i_sensor,3);
        end
    end
end

% %Visual test if projected points  lie on plane
% plot3([0:Vector_InformationLine(1)],[0:Vector_InformationLine(2)], [0:Vector_InformationLine(3)],...
%     '-', 'color', [0.4000    0.8000    1.0000], 'LineWidth', 5)
% hold on
% plot3([0:Vector_DurationLine(1)], [0:0.5:Vector_DurationLine(2)], [0:0.25:Vector_DurationLine(3)], '-',...
%     'color', [1.0000    0.4000    0.4000], 'LineWidth', 5)
% hold on
% fmesh(zplane)
% grid on
% xlabel('x - 150 ms TD')
% ylabel('y - 300 ms TD')
% zlabel('z - 600 ms TD')
% hold on; plot3(Kprime_exp_Original{1}(1:10,1), ...
%     Kprime_exp_Original{1}(1:10,2), ...
%     Kprime_exp_Original{1}(1:10,3),...
%     '.', 'MarkerSize',10, 'MarkerFaceColor','b')
% hold on; plot3(Kprime4AllTD.Vector_2Dplane.perSens.exp{1}(1:10,1), ...
%     Kprime4AllTD.Vector_2Dplane.perSens.exp{1}(1:10,1), ...
%     Kprime4AllTD.Vector_2Dplane.perSens.exp{1}(1:10,1),...
%     '.', 'MarkerSize',10, 'MarkerFaceColor','r')

%% 4) Compute Vector Norm and Vector Angle for each sensor
vector_InformationLine = param.KprimeComparison.coord_EndPoint_InformationLine;
vector_DurationLine = param.KprimeComparison.coord_EndPoint_DurationLine;

for i_win = 1:size(Kprime4AllTD.Coord.perSens.x_exp,2) %3 Time windows
    %Experimental data
    for i_sens = 1:size(Kprime4AllTD.Coord.perSens.x_exp{i_win},2)
        %Compute vector norm (Euclidean length) of the vector
        %connecting origin (0,0,0) and data point (k' for TD: 150,300,600)
        Kprime4AllTD.Norm_2Dplane.perSens.exp{i_win}(i_sens,1) = ...
            norm(Kprime4AllTD.Vector_2Dplane.perSens.exp{i_win}(i_sens,:));
        
        %Compute vector angle to Information line
        Kprime4AllTD.AngleInformationLine_2Dplane.perSens.exp{i_win}(i_sens,1) = ...
            atan2(norm(cross(Kprime4AllTD.Vector_2Dplane.perSens.exp{i_win}(i_sens,:),vector_InformationLine)),...
            dot(Kprime4AllTD.Vector_2Dplane.perSens.exp{i_win}(i_sens,:),vector_InformationLine));
        
        %Compute vector angle to Duration line
        Kprime4AllTD.AngleDurationLine_2Dplane.perSens.exp{i_win}(i_sens,1) = ...
            atan2(norm(cross(Kprime4AllTD.Vector_2Dplane.perSens.exp{i_win}(i_sens,:),vector_DurationLine)),...
            dot(Kprime4AllTD.Vector_2Dplane.perSens.exp{i_win}(i_sens,:),vector_DurationLine));
        
        %Shuffled data
        for i_rep = 1:size(Kprime4AllTD.Coord.perSens.x_null{1},1)
            %Compute vector norm (Euclidean length) of the vector
            %connecting origin (0,0,0) and data point (k' for TD: 150,300,600)
            Kprime4AllTD.Norm_2Dplane.perSens.shuff{i_win}(i_sens,i_rep) = ...
                norm(Kprime4AllTD.Vector_2Dplane.perSens.shuff{i_win}(i_sens,:,i_rep));
            
            %Compute vector angle to Information line
            Kprime4AllTD.AngleInformationLine_2Dplane.perSens.shuff{i_win}(i_sens,i_rep) = ...
                atan2(norm(cross(Kprime4AllTD.Vector_2Dplane.perSens.shuff{i_win}(i_sens,:,i_rep),vector_InformationLine)),...
                dot(Kprime4AllTD.Vector_2Dplane.perSens.shuff{i_win}(i_sens,:,i_rep),vector_InformationLine));
            
            %Compute vector angle to Duration line
            Kprime4AllTD.AngleDurationLine_2Dplane.perSens.shuff{i_win}(i_sens,i_rep) = ...
                atan2(norm(cross(Kprime4AllTD.Vector_2Dplane.perSens.shuff{i_win}(i_sens,:,i_rep),vector_DurationLine)),...
                dot(Kprime4AllTD.Vector_2Dplane.perSens.shuff{i_win}(i_sens,:,i_rep),vector_DurationLine));
        end
    end
end

%% 5) Average vector norm and angle across sensors in specific cluster
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
                    Kprime4AllTD.Norm_2Dplane.(clusterlabel).exp{i_clusterwindow}{i_clusternum} = ...
                        mean(Kprime4AllTD.Norm_2Dplane.perSens.exp{i_clusterwindow}...
                        (ClusterSOI.(param.KprimeComparison.ClusterSOI{i_clustertype})...
                        {i_clusterwindow}{i_clusternum}'));
                    %vector angle to Information line
                    Kprime4AllTD.AngleInformationLine_2Dplane.(clusterlabel).exp{i_clusterwindow}{i_clusternum} = ...
                        mean(Kprime4AllTD.AngleInformationLine_2Dplane.perSens.exp{i_clusterwindow}...
                        (ClusterSOI.(param.KprimeComparison.ClusterSOI{i_clustertype})...
                        {i_clusterwindow}{i_clusternum}'));
                    %vector angle to Duration line
                    Kprime4AllTD.AngleDurationLine_2Dplane.(clusterlabel).exp{i_clusterwindow}{i_clusternum} = ...
                        mean(Kprime4AllTD.AngleDurationLine_2Dplane.perSens.exp{i_clusterwindow}...
                        (ClusterSOI.(param.KprimeComparison.ClusterSOI{i_clustertype})...
                        {i_clusterwindow}{i_clusternum}'));
                    
                    %Shuffled data
                    %Vector norm
                    Kprime4AllTD.Norm_2Dplane.(clusterlabel).shuff{i_clusterwindow}{i_clusternum} = ...
                        mean(Kprime4AllTD.Norm_2Dplane.perSens.shuff{i_clusterwindow}...
                        (ClusterSOI.(param.KprimeComparison.ClusterSOI{i_clustertype})...
                        {i_clusterwindow}{i_clusternum}',:));
                    %vector angle to Information line
                    Kprime4AllTD.AngleInformationLine_2Dplane.(clusterlabel).shuff{i_clusterwindow}{i_clusternum} = ...
                        mean(Kprime4AllTD.AngleInformationLine_2Dplane.perSens.shuff{i_clusterwindow}...
                        (ClusterSOI.(param.KprimeComparison.ClusterSOI{i_clustertype})...
                        {i_clusterwindow}{i_clusternum}',:));
                    %vector angle to Duration line
                    Kprime4AllTD.AngleDurationLine_2Dplane.(clusterlabel).shuff{i_clusterwindow}{i_clusternum} = ...
                        mean(Kprime4AllTD.AngleDurationLine_2Dplane.perSens.shuff{i_clusterwindow}...
                        (ClusterSOI.(param.KprimeComparison.ClusterSOI{i_clustertype})...
                        {i_clusterwindow}{i_clusternum}',:));
                end
            end
        end
    end
end

%% 6) Compute nonparametrical p-values for Norm and Angle for experimental data
for i_win = 1:size(Kprime4AllTD.Coord.perSens.x_exp,2) %3 Time windows
    for i_sens = 1:size(Kprime4AllTD.Coord.perSens.x_exp{i_win},2)
        %per sensor
        %Vector norm (p-val determined by shuff > exp)
        Kprime4AllTD.Norm_2Dplane.perSens.pval_expvsshuff{i_win}(i_sens,1) = ...
            sum(Kprime4AllTD.Norm_2Dplane.perSens.shuff{i_win}(i_sens,:) >= ...
            (Kprime4AllTD.Norm_2Dplane.perSens.exp{i_win}(i_sens))) ...
            / size(Kprime4AllTD.Coord.perSens.x_null{1},1);
        %Vector angle to Information line (p-val determined by shuff < exp)
        Kprime4AllTD.AngleInformationLine_2Dplane.perSens.pval_expvsshuff{i_win}(i_sens,1) = ...
            sum(Kprime4AllTD.AngleInformationLine_2Dplane.perSens.shuff{i_win}(i_sens,:) <= ...
            (Kprime4AllTD.AngleInformationLine_2Dplane.perSens.exp{i_win}(i_sens))) ...
            / size(Kprime4AllTD.Coord.perSens.x_null{1},1);
        %Vector angle to Duration line (p-val determined by shuff < exp)
        Kprime4AllTD.AngleDurationLine_2Dplane.perSens.pval_expvsshuff{i_win}(i_sens,1) = ...
            sum(Kprime4AllTD.AngleDurationLine_2Dplane.perSens.shuff{i_win}(i_sens,:) <= ...
            (Kprime4AllTD.AngleDurationLine_2Dplane.perSens.exp{i_win}(i_sens))) ...
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
                    Kprime4AllTD.Norm_2Dplane.(clusterlabel).pval_expvsshuff{i_clusterwindow}{i_clusternum} = ...
                        sum(Kprime4AllTD.Norm_2Dplane.(clusterlabel).shuff{i_clusterwindow}{i_clusternum} >= ...
                        (Kprime4AllTD.Norm_2Dplane.(clusterlabel).exp{i_clusterwindow}{i_clusternum})) ...
                        / size(Kprime4AllTD.Coord.perSens.x_null{1},1);
                    %vector angle to Information line (p-val determined by shuff < exp)
                    Kprime4AllTD.AngleInformationLine_2Dplane.(clusterlabel).pval_expvsshuff{i_clusterwindow}{i_clusternum} = ...
                        sum(Kprime4AllTD.AngleInformationLine_2Dplane.(clusterlabel).shuff{i_clusterwindow}{i_clusternum} <= ...
                        (Kprime4AllTD.AngleInformationLine_2Dplane.(clusterlabel).exp{i_clusterwindow}{i_clusternum})) ...
                        / size(Kprime4AllTD.Coord.perSens.x_null{1},1);
                    %vector angle to Duration line (p-val determined by shuff < exp)
                    Kprime4AllTD.AngleDurationLine_2Dplane.(clusterlabel).pval_expvsshuff{i_clusterwindow}{i_clusternum} = ...
                        sum(Kprime4AllTD.AngleDurationLine_2Dplane.(clusterlabel).shuff{i_clusterwindow}{i_clusternum} <= ...
                        (Kprime4AllTD.AngleDurationLine_2Dplane.(clusterlabel).exp{i_clusterwindow}{i_clusternum})) ...
                        / size(Kprime4AllTD.Coord.perSens.x_null{1},1);
                end
            end
        end
    end
end

% %% 7) Save output data
savefile = [path_outputdata 'Group_KprimeNormAngle_2Dplane.mat'];
save(savefile, ...
    'Kprime4AllTD', ...
    '-v7.3');

%% 8) Plot data on 2D plane
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

for i_clustertype = 1:length(param.KprimeComparison.ClusterSOI)
    curr_clustertyp = param.KprimeComparison.ClusterSOI{i_clustertype};
    label_curr_clustertyp = ['per' curr_clustertyp];
    disp(['Plotting 2D representations for Kprime data computed for cluster: ' curr_clustertyp])
    
    for i_win = 1:length(ClusterSOI.(curr_clustertyp))
        
        %Determine TW parameters for text plotting
        w1 = num2str(1000 * windows(i_win, 1) / samplingFreq , 3 );
        w2 = num2str(1000 * windows(i_win, 2) / samplingFreq , 3 );
        win_text = [w1 '-' w2 'ms'];
        
        for i_cluster = 1:length(ClusterSOI.(param.KprimeComparison.ClusterSOI{i_clustertype}){i_win})
            
            indices_sensorcluster = ClusterSOI.(param.KprimeComparison.ClusterSOI{i_clustertype})...
                {i_win}{i_cluster}';
            
            if ~isempty(indices_sensorcluster)
                
                %Plot
                h = figure;
                set(gcf,'units','normalized','outerposition',[0 0 1 1])
                set(gcf,'renderer','painters');
                set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
                
                xlim([0 8])
                ylim([0 8])
                view(0,90)
                caxis([0 200])
                grid on
                set(gca,'XTick', [0 2 4 6 8])
                set(gca,'YTick', [0 2 4 6 8])
                set(gca,'ZTick', [0 2 4 6 8])
                
                %Plot orientation lines
                hold on
                plot3([0:Vector_InformationLine(1)],[0:Vector_InformationLine(2)], [0:Vector_InformationLine(3)], ...
                    '-', 'color', [0.4000    0.8000    1.0000], 'LineWidth', 5)
                hold on
                plot3([0:Vector_DurationLine(1)], [0:0.5:Vector_DurationLine(2)], [0:0.25:Vector_DurationLine(3)], ...
                    '-', 'color', [1.0000    0.4000    0.4000], 'LineWidth', 5)
                hold on
                
                %Plot center of mass (across sensors) for shuffled Kprime values for each repetition
                scatter3(...
                    mean(Kprime4AllTD.Vector_2Dplane.perSens.shuff{i_win}(indices_sensorcluster,1,:),1), ...
                    mean(Kprime4AllTD.Vector_2Dplane.perSens.shuff{i_win}(indices_sensorcluster,2,:),1), ...
                    mean(Kprime4AllTD.Vector_2Dplane.perSens.shuff{i_win}(indices_sensorcluster,3,:),1), ...
                    45,  [0.2 0.2 0.2], 'filled', ...
                    'd', 'MarkerFaceAlpha', 0.5);
                %Plot border around shuffled Kprime value distribution
                Xcoord_allShuff = [];
                Ycoord_allShuff = [];
                Zcoord_allShuff = [];
                for i_repetition = 1:length(Kprime4AllTD.Coord.perSens.x_null{i_win})
                    Xcoord_allShuff = [Xcoord_allShuff; Kprime4AllTD.Vector_2Dplane.perSens.shuff{i_win}(indices_sensorcluster,1,i_repetition)];
                    Ycoord_allShuff = [Ycoord_allShuff; Kprime4AllTD.Vector_2Dplane.perSens.shuff{i_win}(indices_sensorcluster,2,i_repetition)];
                    Zcoord_allShuff = [Zcoord_allShuff; Kprime4AllTD.Vector_2Dplane.perSens.shuff{i_win}(indices_sensorcluster,3,i_repetition)];
                end
                border_pointcloudSHUFFclusterSOI = ...
                    boundary([Xcoord_allShuff, Ycoord_allShuff, Zcoord_allShuff]);
                hold on;
                trisurf(border_pointcloudSHUFFclusterSOI, ...
                    Xcoord_allShuff, ...
                    Ycoord_allShuff, ...
                    Zcoord_allShuff, ...
                    'FaceColor', [0.2 0.2 0.2], 'FaceAlpha', 0.1, 'EdgeAlpha', 0.01);
                
                %Plot experimental Kprime values for respective selected sensors
                scatter3(...
                    Kprime4AllTD.Vector_2Dplane.perSens.exp{i_win}(indices_sensorcluster,1), ...
                    Kprime4AllTD.Vector_2Dplane.perSens.exp{i_win}(indices_sensorcluster,2), ...
                    Kprime4AllTD.Vector_2Dplane.perSens.exp{i_win}(indices_sensorcluster,3), ...
                    45,  [0.109 0.596, 0.086], 'filled', ...
                    'MarkerFaceAlpha', 0.5);
                %Plot center of mass across sensors for experimental Kprime distribution
                scatter3(...
                    mean(Kprime4AllTD.Vector_2Dplane.perSens.exp{i_win}(indices_sensorcluster,1),1), ...
                    mean(Kprime4AllTD.Vector_2Dplane.perSens.exp{i_win}(indices_sensorcluster,2),1), ...
                    mean(Kprime4AllTD.Vector_2Dplane.perSens.exp{i_win}(indices_sensorcluster,3),1), ...
                    250,  [0.109 0.596, 0.086], 'filled', ...
                    'd', 'MarkerFaceAlpha', 1);
                %Plot border around around experimental Kprime value
                border_pointcloudEXPclusterSOI = ...
                    boundary([Kprime4AllTD.Vector_2Dplane.perSens.exp{i_win}(indices_sensorcluster,1), ...
                    Kprime4AllTD.Vector_2Dplane.perSens.exp{i_win}(indices_sensorcluster,2), ...
                    Kprime4AllTD.Vector_2Dplane.perSens.exp{i_win}(indices_sensorcluster,3)]);
                hold on;
                trisurf(border_pointcloudEXPclusterSOI, ...
                    Kprime4AllTD.Vector_2Dplane.perSens.exp{i_win}(indices_sensorcluster,1), ...
                    Kprime4AllTD.Vector_2Dplane.perSens.exp{i_win}(indices_sensorcluster,2), ...
                    Kprime4AllTD.Vector_2Dplane.perSens.exp{i_win}(indices_sensorcluster,3), ...
                    'FaceColor', [0.109 0.596, 0.086], 'FaceAlpha', 0.1, 'EdgeAlpha', 0.01)
                
                title({['Kprime values (exp vs. shuff data) projected on 2Dplane computed for ' ...
                    curr_clustertyp '(Cluster: ' num2str(i_cluster) '); TW: ' win_text],...
                    ['Norm: p = ' num2str(Kprime4AllTD.Norm_2Dplane.(label_curr_clustertyp).pval_expvsshuff{i_win}{i_cluster}), ...
                    '; Angle (Information Line): p = ' num2str(Kprime4AllTD.AngleInformationLine_2Dplane.(label_curr_clustertyp).pval_expvsshuff{i_win}{i_cluster}), ...
                    '; Angle (Duration Line): p = ' num2str(Kprime4AllTD.AngleDurationLine_2Dplane.(label_curr_clustertyp).pval_expvsshuff{i_win}{i_cluster})]}, ...
                    'Interpreter','none')
            
                if param.plot.save   == 1
                    path_figs = [paths_NASTD_MEG.Current_outputfig 'ExpKprime/TDcomparison/2Drepresentation/'];
                    mkdir(path_figs)

                    filename     = ['Group_2DKprime_' curr_clustertyp '_Cluster' num2str(i_cluster) '_' win_text '.png'];
                    figfile      = [path_figs filename];
                    saveas(gcf, [figfile], 'png'); %save png version
                    delete(h);
                end  
            end
        end
    end
end

end
