clear

n = 10;
beta = 1.5;
H = (beta-1)/2;
k = 50;
val = 1;
del = 1e-3;
sigma2 = 1;
nSeries = 1e4;

% beta = 1.01
% val         = [0      .5      1        1.5]
% B_pred_mean = [.0014  .4984   1.0129   1.5077]

% beta = 1.5
% val         = [0       1    ]
% B_pred_mean = [-.0011  .9916]

ind = 0;
while ind < nSeries
    [B x] = synthfbmcircul2(2^n+1, H, sigma2);
    
    j = k+1;
    found = 0;
    while j <= length(B)
        if abs( B(j) - val ) < del
            ind = ind + 1;    
            Bs(ind,:) = B(j-(k-1):j);
            x_pred(ind) = fgn_pred(x(j:-1:j-(k-1)), H, sigma2);
            B_pred(ind) = x_pred(ind) + Bs(ind,end);
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
    plot(1:50, Bs(h,:), 'b-');
    plot(51, B_pred(h), 'r.');
end
title(['beta = ' num2str(beta)])

B_pred_mean = mean(B_pred)