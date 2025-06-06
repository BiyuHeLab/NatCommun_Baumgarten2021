function resize_fig(h, width, height)
% resize_fig(h, width, height)
% 
% Change size of figure with handle h.
% width and height are fractions relative to the total screen size.

coord = get(0, 'ScreenSize');
set(h, 'Position', [1 1 width * coord(3), height * coord(4)])

end