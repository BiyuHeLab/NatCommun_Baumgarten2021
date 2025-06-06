function z = rho2z(rho, n)

% Pearson r2z
z = 0.5 * (log(1+rho) - log(1-rho));

% adjustment for Spearman's rho
%z = z * sqrt( (n-3) / 1.06 );