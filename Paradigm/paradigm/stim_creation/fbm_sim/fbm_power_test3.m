clear

N  = 2^10+1;
Hs = [.01, .1:.1:.9, .99];

nreps = 1e2;

use_fGn = 1;

for i_H = 1:length(Hs)

i_H

H = Hs(i_H);

for i_rep = 1:nreps
    
[B x]  = synthfbmcircul2(N, H, 1);

if use_fGn
    % power spectrum estimation on fGn
    beta_PS          = powerspectrum_fft(x,1,0);
    H_PS(i_rep, i_H) = (beta_PS+2-1)/2;

    % DFA on fGn
    alpha             = DFA_copy(x, 0);
    H_DFA(i_rep, i_H) = alpha;

    % wavelet on fGn
    H_wav(i_rep, i_H) = wavelet_est(x,1);

else
    % power spectrum estimation on fBm
    beta_PS          = powerspectrum_fft(B,1,0);
    H_PS(i_rep, i_H) = (beta_PS-1)/2;

    % DFA on fBm
    alpha             = DFA_copy(B, 0);
    H_DFA(i_rep, i_H) = alpha - 1;

    % wavelet on fGn
    H_wav(i_rep, i_H) = wavelet_est(B,0);
end


end

end


%% plot

h = figd(20);
hold on;

set(h, 'PaperPositionMode', 'auto')
set(h, 'Position', [0 0 1500 900])

errorbar(Hs, mean(H_PS), std(H_PS), 'b-');
errorbar(Hs, mean(H_DFA), std(H_DFA), 'r-'); 
errorbar(Hs, mean(H_wav), std(H_wav), 'g-');

plot(Hs, Hs, 'k-');

xlabel('generating H for synthfbmcircul');
ylabel({'estimated H', ['(mean for ' num2str(nreps) ' samples, length = ' num2str(N) ')']});

if use_fGn, mstr = 'fGn'; else mstr = 'fBm'; end
title(['Estimating H using the ' mstr]);

legend('PSD', ...
       'DFA', ...
       'wavelet', ...
       'actual', ...
       'Location', 'NorthWest');


columnname   = {'<html><font size=+2>mean bias (%)', '<html><font size=+2>mean st dev (%)'};
columnformat = {'char', 'numeric', 'numeric'}; 
rowname      = {'<html><font size=+2>PSD', ...
                '<html><font size=+2>DFA', ...
                '<html><font size=+2>wavelet'};

dat =  {mean(abs(mean(H_PS)-Hs))*100, mean(std(H_PS))*100; ...
        mean(abs(mean(H_DFA)-Hs))*100, mean(std(H_DFA))*100; ...
        mean(abs(mean(H_wav)-Hs))*100, mean(std(H_wav))*100 };

t = uitable('Units', 'normalized', ...
            'Data', dat, ... 
            'ColumnName', columnname, ...
            'ColumnFormat', columnformat, ...
            'RowName', rowname, ...
            'FontSize', 15, ...
            'Position', [.45 .14 .45 .14]); 