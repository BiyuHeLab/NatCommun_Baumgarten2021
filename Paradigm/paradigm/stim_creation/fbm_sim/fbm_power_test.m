clear

N     = 2^10+1;
betas = linspace(1.1,2.9,5);

nreps = 1e2;

for i_b = 1:length(betas)

i_b
    
for i_rep = 1:nreps
    
beta  = betas(i_b);
H_fBm = (beta-1)/2;

[B x]  = synthfbmcircul2(N, H_fBm, 1);
beta_PS(i_rep) = powerspectrum_fft(B,1,0);

alpha_B = DFA_copy(B, 0);
% alpha_x = DFA_copy(x);

beta_DFA_B(i_rep) = 2*alpha_B-1;
% beta_DFA_x(i_rep) = 2*alpha_x-1;

end

beta_PS = (beta_PS-1)/2;
beta_DFA_B = (beta_DFA_B-1)/2;

mbeta_PS(i_b) = mean(beta_PS);
mbeta_DFA_B(i_b) = mean(beta_DFA_B);
% mbeta_DFA_x(i_b) = mean(beta_DFA_x);

sbeta_PS(i_b) = std(beta_PS);
sbeta_DFA_B(i_b) = std(beta_DFA_B);
% sbeta_DFA_x(i_b) = std(beta_DFA_x);

% sbeta_PS(i_b) = std(beta_PS) / sqrt(nreps);
% sbeta_DFA_B(i_b) = std(beta_DFA_B) / sqrt(nreps);
% sbeta_DFA_x(i_b) = std(beta_DFA_x) / sqrt(nreps);

end

betas = (betas-1)/2;

fs = 15;
lw = 2;
ms = 10;

figure; hold on;
errorbar(betas, mbeta_PS,sbeta_PS,'bo-','LineWidth',lw,'MarkerSize',ms);
errorbar(betas, mbeta_DFA_B,sbeta_DFA_B,'r^-','LineWidth',lw,'MarkerSize',ms);
% errorbar(betas, mbeta_DFA_x,sbeta_DFA_x,'gV-','LineWidth',lw,'MarkerSize',ms);
plot(betas,betas,'k-','LineWidth',lw);

xlabel('generating H_{fBm} for synthfbmcircul','FontSize',fs)
% ylabel(['estimated H_{fBm} (mean for ' num2str(nreps) ' samples, length = ' num2str(N) ')'],'FontSize',fs)
ylabel({'estimated H_{fBm}', ['(mean for ' num2str(nreps) ' samples, length = ' num2str(N) ')']},'FontSize',fs)

legend('PSD','DFA','actual','Location','NorthWest');
set(gca,'FontSize',fs);