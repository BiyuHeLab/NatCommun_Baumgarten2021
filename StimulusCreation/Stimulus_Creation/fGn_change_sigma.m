function [s_adj fGn_adj fGn] = fGn_change_sigma(s, beta, sigma)
% [s_adj fGn_adj fGn] = fGn_change_sigma(s, beta, sigma)
% 
% Given an fGn/fBm series s with parameter beta, return the same series
% with an adjustment such that the standard deviation of its fGn now equals
% sigma

if beta < 1
    fGn = s;
else
    fGn = s(2:end) - s(1:end-1);
end

fGn_adj = fGn * (sigma / std(fGn) );

if beta < 1
    s_adj = fGn_adj;
else
%     s_adj = [0 cumsum(fGn_adj)] + s(1);
    s_adj = [0 cumsum(fGn_adj)];

end


end

