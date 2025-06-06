function save_fig(h, figfile, resize, type, format)
% save_fig(h, figfile, resize, type, format)
% 
% Save Matlab figure to image file.
% 
% inputs
% ------
% h       - handle for relevant figure. default = gcf
% figfile - filename
%           include entire path in figfile to specify a particular directory
%           otherwise, file is saved to current directory
% resize  - if 1, figure is resized to max screen size before saving
%           default = 0
% type    - method of saving the image
%           1 -> saveas(h, figfile, format);
%           2 -> img = getframe(h); imwrite(img.cdata, figfile, format);
%           default = 1
% format  - image format. default = 'png'


%%

if ~exist('h', 'var') || isempty(h)
    h = gcf;
end

if ~exist('resize', 'var') || isempty(resize)
    resize = 0;
end

if ~exist('type', 'var') || isempty(type)
    type = 1;
end

if ~exist('format', 'var') || isempty(format)
    format = 'png';
end



%%
if resize
    screen_coord = get(0,'ScreenSize') - [0 0 5 5];
    set(h, 'Position', screen_coord);
end


switch type
    case 1
        saveas(h, figfile, format);
    
    case 2
        img = getframe(h);
        imwrite(img.cdata, figfile, format);
        
end