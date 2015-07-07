function G = G_Hat_d(X)

% we use var(X,1,2) instead of var(X,0,2) to obtain the true sample estimator, which is slightly biased
% see "Risk and Asset Allocation"- Springer (2005), by A. Meucci

G = var(X,1,2)+ mean(X,2).^2 - mean(X,2);
