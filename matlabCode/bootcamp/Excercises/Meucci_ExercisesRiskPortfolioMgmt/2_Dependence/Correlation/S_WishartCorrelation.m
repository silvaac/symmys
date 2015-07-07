% this script computes the correlation of the first diagonal and
% off-diagonal elements of a 2x2 Wishart distribution as a function of the inputs
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clear;  close all;  clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=[1 1];
nu=15;

rhos=[-.99 : .01 : .99];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analytical
corrs=sqrt(2)*rhos./sqrt(1+rhos.^2);

% "brute force" check
for uu=1:length(rhos)
    rho=rhos(uu);

    Sigma=diag(s)*[1 rho;rho 1]*diag(s);

    % compute expected values of W_xx and W_xy, see (2.227) in "Risk and Asset Allocation - Springer
    E_xx=nu*Sigma(1,1);
    E_xy=nu*Sigma(1,2);

    % compute covariance matrix of W_xx and W_xy, see (2.228) in "Risk and Asset Allocation - Springer
    m=1;n=1;p=1;q=1;
    var_Wxx=nu*(Sigma(m,p)*Sigma(n,q)+Sigma(m,q)*Sigma(n,p));
    m=1;n=2;p=1;q=2;
    var_Wxy=nu*(Sigma(m,p)*Sigma(n,q)+Sigma(m,q)*Sigma(n,p));
    m=1;n=1;p=1;q=2;
    cov_Wxx_Wxy=nu*(Sigma(m,p)*Sigma(n,q)+Sigma(m,q)*Sigma(n,p));

    S_xx_xy=[var_Wxx cov_Wxx_Wxy
        cov_Wxx_Wxy var_Wxy];

    % compute covariance of X_1 and X_2
    S=diag(1./[sqrt(var_Wxx) sqrt(var_Wxy)]) *  S_xx_xy   * diag(1./[sqrt(var_Wxx) sqrt(var_Wxy)]);

    % correlation = covariance
    corrs2(uu)=S(1,2);
end


figure
plot(rhos,corrs,'.')
hold on
plot(rhos,corrs2,'r')
xlabel('input \rho')
ylabel('Wishart correlation')
