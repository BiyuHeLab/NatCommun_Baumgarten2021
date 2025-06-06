function sensor_topo(sensor_labels, data)

if ~exist('data','var') || isempty(data)
    data = zeros(271,1);
end

topo_data = make_ft_struct(data, 'timelock');


cfg = [];
cfg.layout    = 'CTF275.lay';         
cfg.comment   = 'no';

cfg.highlight        = 'on';
cfg.highlightchannel = sensor_labels;
cfg.highlightsymbol  = 'o';
cfg.highlightcolor   = [0 0 1];
cfg.highlightsize    = 9;

% cfg.marker = 'labels';  
% cfg.markerfontsize = 5;

figure;
ft_topoplotER(cfg, topo_data)