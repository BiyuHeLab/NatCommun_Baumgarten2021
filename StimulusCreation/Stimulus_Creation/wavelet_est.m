function H_est = wavelet_est(data, isfGn)
% H_est = wavelet_est(data, isfGn)
%
% - Estimate the Hurst exponent of data series using wavelet estimation
% - Set isfGn = 1 if data is fGn, 0 if data is fBm (defaults to 1)

if ~exist('isfGn','var') || isempty(isfGn)
    isfGn = 1;
end

% number of vanishing moments (of the Daubechies wavelet).
N = 1;

% wavelet normalisation convention:  
% 1 = L1 [1/scale]       (exponent called  zeta_q),     Zeta_q Diagram,
% 2 = L2 [1/sqrt(scale)] (exponent called alpha_q), Logscale_q Diagram
norm = 2;

% sigtype:  0:  Gaussien, bias et variance of the S_q calculated (in wtspecq_statlot)               (G in title)
% 1:  Non-G with finite variance, bias et variance of the S_q estimated (wtspecq_statlot) (NG in title) 
% 2:  Alpha-stable  - not implemented  (alpha S in title).
sigtype = 0;

% the powers used for S_q(j)  => can be a vector
q = 2; %1:5;

% j1:    the lower limit of the scales chosen,  1<= j1 <=scalemax-1, can be a vector
% j2:    the upper limit of the octaves chosen, 2<= j2 <=scalemax,   can be a vector
j1 = 1; % min(v);
j2 = floor(log2(length(data))); % max(v);

printout = 0;

% 1 --> wavelet increment method (missing a subfunction for this tho)
% 0 --> wavelet method
if isfGn
    increment = 0;
else
    increment = 1;
end

slope = MDestimate3(data, N, norm, sigtype, q, j1, j2, printout, increment);

if isfGn
    offset = 1;
else
    offset = 0;
end

if length(q) > 1
    Hests = q.^-1 .* slope + offset;
    H_est = mean(Hests);
else
    H_est = q^-1 * slope + offset;
end

end