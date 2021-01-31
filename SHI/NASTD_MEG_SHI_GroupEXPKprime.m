function NASTD_MEG_SHI_GroupEXPKprime...
    (subs, tonedur_text, ...
    plot_figs, save_figs, ...
    paths_NASTD_MEG)

%Aim: Average k-prime values across subjects and plot resulting grand 
%average as topoplot

%% 0) Specify vars, paths, and setup fieldtrip
addpath(genpath(paths_NASTD_MEG.ScriptsDir));

%Tone duration condition for data load-in
if strcmp(tonedur_text,'0.15')
    tonedur_title = '0.15sTD';
elseif  strcmp(tonedur_text,'0.3')
    tonedur_title = '0.3sTD';
elseif strcmp(tonedur_text,'0.6')
    tonedur_title = '0.6sTD';
end

path_outputdata = [paths_NASTD_MEG.Current_outputdata 'ExpK/' tonedur_title '/Group/'];
mkdir(path_outputdata)

%% 1) Define analysis parameters
%Define windows
win_size    = 30; 
win_overlap = 0; 

samplingFreq = 600;
tonedur_inSecs  = str2num(tonedur_text);
nSamplesPerTone = tonedur_inSecs * samplingFreq;

windows = [1 win_size];
while windows(end,end) < nSamplesPerTone
    windows = [windows; windows(end,:) + (win_size - win_overlap)];
end

if windows(end,end) > nSamplesPerTone
    windows(end,:) = [];
end

windows_inms = (windows / 600) * 1000;
nWindows = size(windows,1);

%% 2) Load subject data and average kprime values across subjects
%Store all subs in common struct
for i_sub = 1:length(subs)
    sub = subs{i_sub};
    path_inputdata = [paths_NASTD_MEG.Current_outputdata 'ExpK/' tonedur_title '/' sub '/'];
    load([path_inputdata sub '_ExpKprime_' tonedur_title '.mat']);
    %Load in combined exp Kvals sets for each subject
    for i_win = 1:nWindows
        %Create summary matrix with k-prime values (row = subject,column = sensor)
        Kprime_AllSub{i_win}(i_sub, :) = Kprime_data.avgFolds.Kprime.Avg{i_win};
    end
    clear Kprime_data
end

%% 3) Change fomat of across-subject k-prime average to FT-structure for later plotting
%find zlimits for constant plotting across time windows
dv_absmax = 0;
dv_absmin = 16;
for i_win = 1:nWindows
    Kprime_Gavg = mean(Kprime_AllSub{i_win});
    %avg across subjects yields 1x272 k-prime values for current model and window
    
    % find zlimits for constant plotting across time windows
    dv_win_absmax = max(abs(Kprime_Gavg));
    dv_win_absmin = min(abs(Kprime_Gavg));
    
    if dv_win_absmax > dv_absmax
        dv_absmax = dv_win_absmax;
    end
    if dv_win_absmin < dv_absmin
        dv_absmin = dv_win_absmin;
    end
    dv_absmax_plot = ceil(dv_absmax);
    dv_absmin_plot = floor(dv_absmin);
    
    
    load([paths_NASTD_MEG.ScriptsDir 'MEG_sensor_setup_272/label272.mat'] ); %file with CTF sensor labels for 272 sensors
    
    data_topoplot{i_win}.dimord     = 'chan_time';
    data_topoplot{i_win}.label      = label;
    data_topoplot{i_win}.time       = 0;
    data_topoplot{i_win}.avg        = Kprime_Gavg';
end

%% 4) Plot across-subject average k-prime values in topoplot for each time window
if plot_figs == 1
    
    %Determine analysis option for title
    model_comp_title = 'Group-level Kprime per sensor';
    
    %Plot k-prime across-subject average for each time window
    for i_win = 1:nWindows
        %determine start & end points for respective window
        w1 = num2str(1000 * windows(i_win, 1) / samplingFreq , 3 );
        w2 = num2str(1000 * windows(i_win, 2) / samplingFreq , 3 );
        win_title = ['TW = [' w1 ' ms - ' w2 ' ms]'];
        
        cfg = [];
        cfg.layout    = 'CTF275.lay';
        cfg.comment   = 'no';
        cfg.colorbar  = 'yes';
        cfg.zlim = [dv_absmin_plot, dv_absmax_plot]; % 'maxabs';
        
        h = figure;
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        ft_topoplotER(cfg, data_topoplot{i_win});
        
        title({[model_comp_title], ['TD = ' tonedur_text '0ms; ERF ' win_title]})
        
        if save_figs    == 1
            path_figs = [paths_NASTD_MEG.Current_outputfig 'ExpKprime/' tonedur_title '/Group/'];
            mkdir(path_figs)
            
            filename     = ['Group_Kprime_' tonedur_title '0_TW' w1 '-' w2 '.png'];
            figfile      = [path_figs filename];
            saveas(gcf, [figfile], 'png'); %save png version
            delete(h);
        end      
    end
end

end
