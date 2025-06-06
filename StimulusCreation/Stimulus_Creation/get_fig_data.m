function [xdata, ydata, zdata] = get_fig_data(h)
% [xdata, ydata, zdata] = get_fig_data(h)
%
% Retrieve the data plotted in Matlab figure h.
%
% via
% http://www.mathworks.com/matlabcentral/answers/100687-how-do-i-extract-data-from-matlab-figures

axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes

objTypes = get(dataObjs, 'Type');  %type of low-level graphics object

xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs, 'YData');
zdata = get(dataObjs, 'ZData');