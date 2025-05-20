H = 1e-2;
d = H - .5;
N = 2^10;

Hs = [.01 .1:.1:.9 .99];

for h=1:length(Hs)
    h
for i=1:1000
[B x]   = synthfbmcircul(N,Hs(h));
s(i,h)  = std(x);
sB(i,h) = std(B);
end
end

figure; hold on;
plot(Hs,mean(s),'bo-')
plot(Hs,mean(sB),'ro-')

N=2^13;
figure;
[B1 x1] = synthfbmcircul(N,.01);
[B2 x2] = synthfbmcircul(N,.5);
x2 = x2 * 2;
B2 = [0 cumsum(x2(1:N-1))];
subplot(211); hold on;
plot(B1,'b-');
plot(B2,'r-');
subplot(212); hold on;
plot(x1,'b-');
plot(x2,'r-');