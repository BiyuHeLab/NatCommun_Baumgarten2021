function [hh hv] = hv_errorbars(x,y,xe,ye,style,lineWidth)

if ~exist('style','var')
    style = 'b-';
end

if ~exist('lineWidth','var')
    lineWidth = 1;
end

for i = 1:length(x)
    hh(i) = plot([x(i)+xe(i) x(i)-xe(i)],[y(i) y(i)],style,'LineWidth',lineWidth);    
    hv(i) = plot([x(i) x(i)],[y(i)+ye(i) y(i)-ye(i)],style,'LineWidth',lineWidth);
end