function NASTD_MEG_Preproc_ICAselection...
    (sub, paths_NASTD_MEG)

%Aim: Select ICA-components containing artifacts

%% 0) Specify vars, paths, and setup fieldtrip
addpath(genpath(paths_NASTD_MEG.ScriptsDir));
addpath(paths_NASTD_MEG.FieldTrip);
ft_defaults

NASTD_MEG_SubInfo

loadfile_ica = [si.path_ica  'data_ica.mat'];

%% 1) Load preprocessed data

disp('loading...')
load(loadfile_ica);
disp('done.')


%% 2) make a topological plot of the independent components 
% save the identity of the artifactual ICs in 
% path_ica/ICA_artifact_selection.m

screen_size = get(0,'ScreenSize');
screen_l = screen_size(3);
screen_h = screen_size(4);

comp1 = 11;

close all
% ICA_fft_subplot(data_ica,4,5,comp1);

cfg           = [];
cfg.component = comp1 : comp1 + 9;        % specify the component(s) that should be plotted
cfg.layout    = 'CTF275.lay'; % specify the layout file that should be used for plotting
cfg.comment   = 'no';

ft_topoplotIC(cfg, data_ica);
set(gcf, 'Position', [1, 1, (screen_l/2)-1, screen_h]);

cfg           = [];
cfg.layout    = 'CTF275.lay';
cfg.viewmode  = 'component';

ft_databrowser(cfg, data_ica);
set(gcf, 'Position', [(screen_l/2), 1, screen_l/2, screen_h]);

ICA_fft_subplot(data_ica,3,4,comp1);

return

%% 3) Beginn manual selection 
% double-check selected components
% this section of code is meant to be run manually as a follow-up to manual
% inspection of the ICA data for artifact selection

% manually selected artifacts to be removed

run([si.path_ica 'selected_artifacts.m']);
artifact_components = ac_all;

for i = 1:length(artifact_components)
    comp_label{i,1} = data_ica.label{artifact_components(i)};
end

cfg           = [];                                           
cfg.component = artifact_components;
cfg.layout    = 'CTF275.lay';         
cfg.comment   = 'no';

ft_topoplotIC(cfg, data_ica);
set(gcf, 'Position', [1, 1, (screen_l/2)-1, screen_h]);

cfg           = [];
cfg.channel   = comp_label;
cfg.layout    = 'CTF275.lay';
cfg.viewmode  = 'component';

ft_databrowser(cfg, data_ica);
set(gcf, 'Position', [(screen_l/2), 1, screen_l/2, screen_h]);

return

