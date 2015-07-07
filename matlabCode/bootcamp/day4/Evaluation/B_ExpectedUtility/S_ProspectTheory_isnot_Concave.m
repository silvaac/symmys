% this script shows that the certainty equivalent for prospect theory 
% (i.e. "erf" utility function) investor is not a concave function of allocation
% see "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clear; clc; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs (normal market, prospect-theory utility)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mu=[1 1]';
sigma=[1 1];
Rho=-.99;

Alpha_1s=[-1:.1:1];
Alpha_2s=[-1:.1:1];

NumSimulations=10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=[1 Rho; Rho 1];
Sigma=diag(sigma)*C*diag(sigma);


CE=zeros(length(Alpha_1s),length(Alpha_2s));
for n=1:length(Alpha_1s)
    Countdown = length(Alpha_1s)-n+1
    for m=1:length(Alpha_2s)
        Alpha=[Alpha_1s(n); Alpha_2s(m)];
        mu=Alpha'*Mu;
        sigma=sqrt(Alpha'*Sigma*Alpha);
        Portfolio=mu+sigma*randn(NumSimulations,1);
        Exp_U=mean(erf(Portfolio));
        CE(n,m)=erfinv(Exp_U);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ALPHA_1s,ALPHA_2s] = meshgrid(Alpha_1s,Alpha_2s);
figure
hf=surf(ALPHA_1s,ALPHA_2s,CE');
grid on


