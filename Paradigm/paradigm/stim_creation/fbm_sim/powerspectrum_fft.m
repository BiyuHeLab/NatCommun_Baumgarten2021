function beta = powerspectrum_fft(timeseries,T,doplot)

if ~exist('doplot','var'), doplot = 1; end

n=length(timeseries);
p=abs(fft(timeseries))/(n/2);

n2_int = floor(n/2);
p=p(1:n2_int).^2;
freq=[0:n2_int-1]/T;

%  semilogy(freq,p)

logf = log10(freq(2:end));
logp = log10(p(2:end));

c  = polyfit(logf, logp, 1);
fc = [min(logf) max(logf)];
pc = c(1)*fc+c(2);

beta = -c(1);


if doplot
    
figure; 
subplot(1,2,1);
plot(timeseries);

subplot(1,2,2);
hold on;
plot(freq,p);

plot(10.^fc, 10.^pc, 'r-');
title(['beta = ' num2str(beta)]);

set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log'); 

end