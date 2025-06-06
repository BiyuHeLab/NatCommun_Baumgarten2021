function smoothed = myft_smooth(sourcespace, unit, fwhm_inmm)
% smoothed = myft_smooth(sourcespace, unit, fwhm_inmm)
%
% sourcespace is a struct holding the data to be smoothed.
% it should have fields
% - pos : Nx3 matrix of spatial location of each point in the mesh
% - avg.pow : Nx1 vector of activity at each position
%
% unit indicates the units of sourcespace.pos. 
% for the typical processing pipeline this will be 'cm'
%
% fwhm_inmm is the FWHM of the smoothing, specified in mm.
% if the "unit" input is 'cm' rather than 'mm' then this function
% automatically converts fwhm_inmm to cm.

%%

s = sourcespace;
s.unit = unit;

if strcmp(unit, 'cm')
    s.fwhm = fwhm_inmm / 10;
    
elseif strcmp(unit, 'mm')
    s.fwhm = fwhm_inmm;
    
else
    disp('unrecognized units!')
    return

end


%% 

% conversion from FWHM to sigma taken from
% http://support.brainvoyager.com/functional-analysis-preparation/27-pre-processing/279-spatial-smoothing-in-preparation.html
sigma = s.fwhm / sqrt(8*log(2));
sigma = sigma * eye(3);


%%

npts = size(s.pos, 1);
nt   = size(s.avg.pow, 2);

pow_smooth = zeros(npts, nt);
for i_pt = 1:npts
    mu = s.pos(i_pt, :);
    w  = mvnpdf(s.pos, mu, sigma);
    w  = w / sum(w);

%     if i_pt==121, keyboard, end
    
    for i_t = 1:nt
        pow_smooth(i_pt, i_t) = sum( w .* s.avg.pow(:, i_t) );
    end
end


smoothed = s;
smoothed.avg.pow = pow_smooth;

    
