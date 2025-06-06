function xd = discretize(x, vals)
% xd = discretize(x, vals)
% 
% for each element v in the vector x, set v to the closest value contained in vals

xd = zeros(size(x));
for j = 1:length(x)
    diff    = abs(x(j) - vals);
    [m ind] = min(diff);
    xd(j)   = vals(ind);
end