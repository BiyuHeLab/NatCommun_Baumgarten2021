N     = 2^13;
beta  = 1.3;
H_fBm = (beta-1)/2;

[B x]  = synthfbmcircul(N, H_fBm);

beta_true = beta
beta_PS = powerspectrum_fft(B,1,1)

alpha_B = DFA_copy(B);
alpha_x = DFA_copy(x);

beta_DFA_B = 2*alpha_B-1
beta_DFA_x = 2*alpha_x-1