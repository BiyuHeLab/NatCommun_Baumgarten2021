%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    synthfbmcircul.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Fractional Brownian motion synthesis with circulant matrix method
% B(t) is synthetized for t in [0,1] and variance (of white gaussian noise) 
% equal to 1.
%
% Inputs:
%          N : number of samples
%          H : Hurst parameter (0 < H < 1)
%     [seed] : states (integers) for random number
%              generator (real sequence and imaginary sequence)
% Outputs:
%          B : H-fBm N-length trace
%          x : H-fGn : x[n] = B[n+1] - B[n]
%          w : white gaussian noise generator (complex)
%
%  >  > >>  Take N = 2^k + 1  << <  <
%
% Example :
%-----------
% N=2^10+1 ;   % N=1025
% H=1/3 ;
% [B,x,w] = synthfbmcircul(N,H) ;

function [B,x,w] = synthfbmcircul(N,H,seed) ;
B = 0; x =0 ; w= 0 ;

if exist('seed','var')
    randn('state',seed)
end

s = 1 ;
n = 0:N-1 ;
tmax = 1 ;
dt = tmax/N ;

% Synthese de la covariance du fGn

r = (dt)^(2*H)*s/2*(abs(n-1).^(2*H) + abs(n+1).^(2*H) - 2*abs(n).^(2*H)) ;

% Insertion dans matrice circulante :

rz = [r fliplr(r(2:N-1))] ;
clear r

% Calcul des racines de la matrice circulante (=Fourier sur rz)
L = real(fft(rz)) ;
vmin=min(min(real(L))') ;
% sum(abs(w) - w)
clear rz
% Generation bruit blanc

wr =  randn(1,2*N-2) ;

wi = randn(1,2*N-2) ;
w = wr + i.*wi ;

clear wr wi
w2 = sqrt(L./(2*N-2)).*w ;


clear L w
%% ATTENTION: ne pas utiliser ifft, qui utilise une normalisation differente

z = fft(w2) ;
clear w2
% Processus increment normalise

x = real(z(1:N)) ;
clear z

% Processus integre: fBm

B = [0 cumsum(x(1:N-1))] ;
