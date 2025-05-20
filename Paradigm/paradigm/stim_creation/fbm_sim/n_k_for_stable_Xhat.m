clear

Hs    = [1e-2 .25 .5 .7 .9 1-1e-2];
N     = 2^10+1;
sigma = 1;

ks     = 1:55;
nreps  = 25;


for i_H = 1:length(Hs)
    i_H
    
    H = Hs(i_H);
    d = H - .5;

    for i_rep = 1:nreps

        [B x] = synthfbmcircul2(N,H,sigma);
        Xt(i_H, i_rep) = x(end);

        for i_k = 1:length(ks)
            Xts = x(end-1 : -1 : end - ks(i_k));
            Xt_est(i_H, i_rep, i_k) = arima_pred(Xts,d);
        end
    end

end

% plot
fs = 15;

figure;
for i_H = 1:length(Hs)
    subplot(3,2,i_H); hold on;
    for i_rep = 1:nreps
        plot(ks, Xt(i_H,i_rep) - squeeze(Xt_est(i_H, i_rep, :)),'b-');
    end
    set(gca,'FontSize',fs);
    
    title(['H = ' num2str(Hs(i_H))],'FontSize',fs);
    ylabel('X_t - X_{t-hat}(d,k)','FontSize',fs);
    xlabel('k','FontSize',fs);

    xlim([min(ks) max(ks)])
end



figure;
for i_H = 1:length(Hs)
    subplot(3,2,i_H); hold on;
    plot(ks, Xt(i_H,1) - squeeze(Xt_est(i_H, 1, :)),'b-');

    set(gca,'FontSize',fs);
    
    title(['H = ' num2str(Hs(i_H))],'FontSize',fs);
    ylabel('X_t - X_{t-hat}(d,k)','FontSize',fs);
    xlabel('k','FontSize',fs);

    xlim([min(ks) max(ks)])
end