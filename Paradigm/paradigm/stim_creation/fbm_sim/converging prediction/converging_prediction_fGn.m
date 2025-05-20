clear

n = 10;
beta = .5;
H = (beta+1)/2;
k = 50;
val = 1.5;
del = 1e-3;
sigma2 = 1;
nSeries = 1e4;

% beta = 1.01
% val         = [0      .5      1       1.5]
% x_pred_mean = [.0011  .2039   .4131   .6189]


ind = 0;
while ind < nSeries
    [B x] = synthfbmcircul2(2^n+1, H, sigma2);
    
    j = k+1;
    found = 0;
    while j <= length(x)
        if abs( x(j) - val ) < del
            ind = ind + 1;            
            xs(ind,:) = x(j-(k-1):j);
            x_pred(ind) = fgn_pred(x(j:-1:j-(k-1)), H, sigma2);
%             x_pred(ind)  = x_arima(ind);
            found = 1;
        end
        
        if found
            j = Inf;
        else
            j = j+1;
        end
    end
end

figure; hold on;
for h=1:100
    plot(1:50, xs(h,:), 'b-');
    plot(51, x_pred(h), 'r.');
end
title(['beta = ' num2str(beta)])

x_pred_mean = mean(x_pred)