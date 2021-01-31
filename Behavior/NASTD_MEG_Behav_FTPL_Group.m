function  NASTD_MEG_Behav_FTPL_Group...
    (subs, tonedur_text, ...
    plot_figs, save_figs, ...
    paths_NASTD_MEG)
%Compute and plot behavioral Final Tone Pitch Likelihood (FTPL) data 
%for group-level (Corresponds to Fig. 2)

%% 1) Load behavioral data and check trial numbers
for i_sub = 1:length(subs)
    
    sub = subs{i_sub}
    NASTD_MEG_SubInfo%load subject info file (var: si)
    
    load([paths_NASTD_MEG.Current_inputdata 'subjects/sfa_expt3_' sub '.mat']);
    
    %Filter behavioral data to exclude trials that don't have MEG recording
    numoftrials = num2str(length(data.trialNum));
    
    dat_fields  = fieldnames(data);
    for i = 1:length(dat_fields)
        if eval(['length(data.' dat_fields{i} ') == 324'])
            eval(['data.' dat_fields{i} ' = data.' dat_fields{i} '(si.good_trials);']);
        end
    end
    
    stim_fields = fieldnames(data.stim);
    for i = 1:length(stim_fields)
        if eval(['length(data.stim.' stim_fields{i} ') == 324'])
            eval(['data.stim.' stim_fields{i} ' = data.stim.' stim_fields{i} '(si.good_trials);']);
        end
    end
    
    stim_fields = fieldnames(stim);
    for i = 1:length(stim_fields)
        if eval(['length(stim.' stim_fields{i} ') == 324'])
            eval(['stim.' stim_fields{i} ' = stim.' stim_fields{i} '(si.good_trials);']);
        end
    end
    
    %% 2) define filters for tone-duration specific trial selection
    %General filter - ensure all behavioral responses were entered properly
    filter_resp = data.resp_beta >= 0 & data.resp_prob >= 0 & data.trialNum > 0;
    filter_rt   = data.r1_RT > 0 & data.r2_RT > 0;
    
    %Tone duration filter
    switch tonedur_text
        case '0.15'
            filter_toneDur = stim.toneDur == 0.15;
        case '0.3'
            filter_toneDur = stim.toneDur == 0.3;
        case '0.6'
            filter_toneDur = stim.toneDur == 0.6;
        otherwise
            filter_toneDur = ones(1,length(data.trialNum));
    end
    
    %Manually defined filters
    filter_manual = ones(1,length(data.trialNum));
    
    %Combined filter
    filter = filter_resp & filter_rt & filter_toneDur & filter_manual;
    
    %% 3) apply filters
    %On MEG data
    dfn = fieldnames(data);
    for j = 1:length(dfn)
        if eval(['length(data.' dfn{j} ') == length(data.trialNum)'])
            eval(['dataf.' dfn{j} ' = data.' dfn{j} '(filter);']);
        end
    end
    
    %On stim subfield (stimulus-wise parameters)
    sfn = fieldnames(stim);
    for j = 1:length(sfn)
        if eval(['length(stim.' sfn{j} ') == length(data.trialNum)'])
            eval(['stimf.' sfn{j} ' = stim.' sfn{j} '(filter);']);
        end
    end
    
    %% 4) Final Tone Probability Rating
    p34 = unique(stimf.logf_final); %possible log freq of final tone (p34)
    predp34ID = [-1 0 1]; %numerical label for predicted final tone (p*34; low/middle/high)
    
    for i_pred = 1:3 %index for predicted final tone
        
        filter_predp34 = stimf.predID == predp34ID(i_pred); %filter to select those trials with specific final tone
        
        for i_p34 = 1:length(p34) %loop across all actual final tone pitches
            
            filter_p34  = stimf.logf_final == p34(i_p34); %filter to select only trials with specific final tone pitch
            
            % get rFTPL rating per presented final tone (p34)
            FTPLrating_perp34(i_sub, i_p34) = ...
                mean(dataf.resp_prob(filter_p34)); 
            SEM_FTPLrating_perp34(i_sub, i_p34) = ...
                std(dataf.resp_prob(filter_p34)) / sqrt(sum(filter_p34));
            
            % get rFTPL rating per predicted final tone (predp34)
            FTPLrating_perpredp34(i_sub, i_pred, i_p34) = ...
                mean(dataf.resp_prob(filter_predp34 & filter_p34));
            SEM_FTPLrating_perpredp34(i_sub, i_pred, i_p34) = ...
                std(dataf.resp_prob(filter_predp34 & filter_p34)) / sqrt(sum(filter_predp34 & filter_p34));
        end
        
        mean_predp34(i_pred) = mean(stimf.logf_pred(filter_predp34));
        
    end
end

%% 5) Average behavioral output variables across subs
FTPLrating_perp34 = mean(FTPLrating_perp34);
SEM_FTPLrating_perp34 = mean(SEM_FTPLrating_perp34);

GAvg_FTPLrating_perpredp34 = squeeze(mean(FTPLrating_perpredp34));
GAvg_SEM_FTPLrating_perpredp34 = squeeze(mean(SEM_FTPLrating_perpredp34));

%% 6) Plotting group-average behavioral results
if plot_figs == 1
    path_fig = [paths_NASTD_MEG.Current_outputfig];
    switch tonedur_text
        case 'all'
            path_fig = [path_fig 'allTD/'];
        case '0.15'
            path_fig = [path_fig '0.15sTD/'];
        case '0.3'
            path_fig = [path_fig '0.3sTD/'];
        case '0.6'
            path_fig = [path_fig '0.6sTD/'];
    end
    
    mkdir(path_fig);
    
    %Plot FTPL rating per p34
    xlab = {'220'  '277'      '349'      '554'      '698'      '880'};
    
    h = figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    
    e = errorbar(p34, FTPLrating_perp34, SEM_FTPLrating_perp34, 'ko-');
    e.LineWidth = 3;
    xlabel('Final tone pitch (Hz)','FontSize',16)
    ylabel('FTPL rating','FontSize',16)
    ylim([1 5])
    xlim([log(200) log(950)])
    set(gca,'XTick',p34)
    set(gca,'XTickLabel', xlab,'FontSize',16)
    set(gca,'YTick',[1:5])
    title({['GAvg (n = ' num2str(length(subs)) ') - FTPLrating per p34 - ToneDur: ' tonedur_text]},'FontSize',20);
    
    if save_figs == 1 %Save Figures?
        filename = [path_fig 'GAvg_FTPLrating_perp34_' tonedur_text 's.png'];
        saveas(gcf, [filename], 'png'); %save png version
        delete(h);
    end
    
    
    %Plot FTPL rating per p34 and p*34
    xlab = {'220'  '277'      '349'      '554'      '698'      '880'};
    
    h = figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    
    etitle = {'p*34 = low', 'p*34 = med', 'p*34 = high'};
    colors = {'bo-' 'co-' 'yo-'};
    colors_line = {'b-' 'c-' 'y-'};
    pred_colors = {'bs', 'cs', 'ys'};
    
    for i_pred = 1:3
        hold on;
        e = errorbar(p34, GAvg_FTPLrating_perpredp34(i_pred,:), GAvg_SEM_FTPLrating_perpredp34(i_pred,:), colors{i_pred}, 'LineWidth', 3, 'MarkerSize', 10);
                
        %Optional: Add individual data points
        for i_sub = 1:length(subs)
            hold on;
            scatter(p34 - ((0.03*i_pred) - 0.06), squeeze(FTPLrating_perpredp34(i_sub,i_pred,:))','o',pred_colors{i_pred},'filled')
%             plot(p34 - ((0.03*i_pred) - 0.06),squeeze(FTPLrating_perpredp34(i_sub,i_pred,:))',colors_line{i_pred},'LineWidth',1) %plot line connecting data points
            
        end
        
        xlabel('Final tone pitch (Hz)','FontSize',16)
        ylabel('FTPL rating','FontSize',16)
        ylim([1 5])
        xlim([log(200) log(950)])
        set(gca,'XTick',p34)
        set(gca,'XTickLabel', xlab,'FontSize',16)
        set(gca,'YTick',[1:5])
        title(['GAvg - FTPLrating per p*34 - ToneDur: ' tonedur_text],'FontSize',20);
        
      
        %         hl = legend(etitle, 'Location', 'NorthEast');
        %         set(hl, 'box', 'off')
        %         set(hl, 'FontSize', 25)
    end
    
    %Add pitch of p*34 on x-axis
    hold on
    for i_pred = 1:3
        plot(mean_predp34(i_pred), 2, pred_colors{i_pred},'LineWidth', 3,'MarkerSize', 20)
        hold on
    end
        
    if save_figs == 1 %Save Figures?
        filename = [path_fig 'GAvg_FTPLrating_perpredp34_' tonedur_text 's.png'];
        saveas(gcf, [filename], 'png'); %save png version
        delete(h);
    end
   
end