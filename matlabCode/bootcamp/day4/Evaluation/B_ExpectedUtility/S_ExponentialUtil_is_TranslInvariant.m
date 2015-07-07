% this script shows that the certainty equivalent for exponential utility functions 
% is translation invariant
% see "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clear; clc; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs (log-normal market, exponential utility)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mu=[0 0]';
sigma=[.000000001 2];
Rho=-.0;
NumSimulations=10000;
Zeta=3;

FixRiskFree_1=.3;
FixRiskFree_2=.8;

Alpha_1s=[0:.05:1];
Alpha_2s=[0:.05:1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=[1 Rho; Rho 1];
Sigma=diag(sigma)*C*diag(sigma);
CompoundedReturns=mvnrnd(Mu,Sigma,NumSimulations);
Prices=exp(CompoundedReturns);

CE=zeros(length(Alpha_1s),length(Alpha_2s));
mm=[]; t=[]; l=[];
for n=1:length(Alpha_1s)
    Countdown = length(Alpha_1s)-n+1
    for m=1:length(Alpha_2s)
        Alpha=[Alpha_1s(n); Alpha_2s(m)];
        Portfolio=Prices*Alpha;
       
        U=-exp(-Portfolio/Zeta);
        Exp_U=mean(U);
        CE(n,m)=-Zeta*log(-Exp_U);
              
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ALPHA_1s,ALPHA_2s] = meshgrid(Alpha_1s,Alpha_2s);
figure
hf=surf(ALPHA_1s,ALPHA_2s,CE');

% plot first plane for fixed risk-free allocation
Zz=[0 max(max(CE)) max(max(CE)) 0];
a2=[Alpha_2s(1) Alpha_2s(1) Alpha_2s(end) Alpha_2s(end)];
a1=FixRiskFree_1+0*a2;
hold on
h=fill3(a1,a2,Zz,Zz);

% plot second plane for fixed risk-free allocation
a1=FixRiskFree_2+0*a2;
hold on
h=fill3(a1,a2,Zz,Zz);

grid on