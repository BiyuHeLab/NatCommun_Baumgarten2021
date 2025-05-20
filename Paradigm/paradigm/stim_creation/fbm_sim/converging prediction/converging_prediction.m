function [series series_pred pred_err fGns fBms] = converging_prediction(beta, sigma2, k, doPlot, nSeries, delta, v_f, progress)
% [series series_pred pred_err] = converging_prediction(beta, sigma2, k, doPlot, nSeries, delta, v_f)
% 
% generate fGn / fBm series that converge to a similar final value, and
% calculate the prediction for the upcoming value
%
% inputs
% ------
% beta     - beta of the fGn / fBm series
% sigma2   - variance of the fGn (default value = 1) 
% k        - # of prior samples on which the prediction for the next sample
%            is made (default = 50)
% doPlot   - set to 1 to plot results, 0 to turn off plot (default = 1)
% nSeries  - # of series to generate (default = 100)
% delta    - tolerance level for final value v_f. i.e. series are only
%            returned if |series(end) - v_f| <= delta (default = 1e-3 * sqrt(sigma2))
% v_f      - value towards which the series converge (default = 1.5 * sqrt(sigma2))
% progress - set to 1 to get a progress report while simulations are
%            ongoing. (default = 0)
%
%
% outputs
% -------
% series      - A matrix of size nSeries x k. Each row is a fGn/fBm of length
%               k for which the final value is within delta of v_f.
% series_pred - A vector of size nSeries x 1. series_pred(i) contains the best 
%               prediction for the value following series(i,k) (found using 
%               function fgn_pred)
% pred_err    - A vector of size nSeries x 1. pred_err(i) contains the
%               prediction error for the k+1'th sample, 
%               i.e. series_pred(i) - series(i,k+1) [series(i,k+1) not
%               included in output]


%% set default input values for unspecified inputs
if ~exist('sigma2','var') || isempty(sigma2)
    sigma2 = 1;
end

if ~exist('k','var') || isempty(k)
    k = 50;
end

if ~exist('doPlot','var') || isempty(doPlot)
    doPlot = 1;
end

if ~exist('nSeries','var') || isempty(nSeries)
    nSeries = 100;
end

if ~exist('delta','var') || isempty(delta)
    delta = 1e-3 * sqrt(sigma2);
end

if ~exist('v_f','var') || isempty(v_f)
    v_f = 1.5 * sqrt(sigma2);
end

if ~exist('progress','var') || isempty(progress)
    progress = 0;
end

%% prepare for search

% ensure synthesized series are long enough to accommodate k
n = 10;
while 2^n + 1 < 5*k
    n = n+1;
end

% find H from beta
if beta > 1
    H = (beta-1)/2;
else
    H = (beta+1)/2;
end


%% search for series that converge to v_f, within tolerance delta

% series counter. 1 <= i_s <= nSeries
i_s = 0;

if progress
    fprintf(1,'\n\npercent done:\n');
end

% while we haven't found the # of required series...
while i_s < nSeries
    
    % synthesize an fGn/fBm
    [B x] = synthfbmcircul2(2^n+1, H, sigma2);

    % select the appropriate series for analysis
    if beta > 1
        s = B;
    else
        s = x;
    end
    
    
    % loop through the series to find a segment of length k for which the
    % k'th value is within delta of v_f
    j = k+1;
    while j+1 <= length(s)
        
        % if we find a qualifying segment of the series...
        if abs( s(j) - v_f ) < delta
            i_s = i_s + 1;
            
            if progress
                percent = 100*i_s/nSeries;
                fprintf(1,'%3.1f%%\n',percent);
            end
            
            
            % store the last k samples of this segment
            series_maxInd = j;
            series_minInd = j - (k-1);
            series(i_s,:)   = s(series_minInd : series_maxInd);
                            
            if beta > 1
                % find the best k'th order prediction for the next fGn sample
                %
                % N.B. that for fBm prediction, the index to the corresponding
                % sample of the fGn is shifted down by 1, since by
                % convention of synthfbmcircul, B[n+1] = B[n] + x[n]

                fGn_maxInd = series_maxInd - 1;
                fGn_minInd = series_minInd - 1;
                xt = x(fGn_maxInd : -1 : fGn_minInd);
                
                x_pred = fgn_pred(xt, H, sigma2);
                series_pred(i_s,1) = x_pred + series(i_s,end);
                
                fBms(i_s,:) = series(i_s,:);
                fGns(i_s,:) = x(fGn_minInd : fGn_maxInd);
            else
                % find the best k'th order prediction for the next fGn sample
                
                fGn_maxInd = series_maxInd;
                fGn_minInd = series_minInd;
                xt = x(fGn_maxInd : -1 : fGn_minInd);
                
                x_pred = fgn_pred(xt, H, sigma2);
                series_pred(i_s,1) = x_pred;
                
                fGns(i_s,:) = series(i_s,:);
                fBms(i_s,:) = B(fGn_minInd+1 : fGn_maxInd+1);
            end
            
            pred_err(i_s,1) = series_pred(i_s,1) - s(j+1);
            
            % exit the search loop
            j = Inf;
        end

        j = j+1;
   end
end


%% plot results

if doPlot
    
    figure; hold on;
    for h = 1:nSeries
        plot(1:k, series(h,:), 'b-');
        plot(k+1, series_pred(h), 'r.');
    end
    title(['beta = ' num2str(beta)])

end