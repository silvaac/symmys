% this script shows that the certainty equivalent for a power utility functions 
% is homogeneous
% see "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clear; clc; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs (log-normal market, power utility)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mu=[1 2]';
sigma=[2 5];
Rho=-.5;
NumSimulations=10000;
Gamma=.03;

Alpha_1s=[0:.05:1];
Alpha_2s=[0:.05:1];

Line=2/3; % radial line

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=[1 Rho; Rho 1];
Sigma=diag(sigma)*C*diag(sigma);
CompoundedReturns=mvnrnd(Mu,Sigma,NumSimulations);
Prices=exp(CompoundedReturns);

% compute CE
CE=zeros(length(Alpha_1s),length(Alpha_2s));
for n=1:length(Alpha_1s)
    Countdown = length(Alpha_1s)-n+1
    for m=1:length(Alpha_2s)
        Alpha=[Alpha_1s(n); Alpha_2s(m)];
        Portfolio=Prices*Alpha;
       
        U=Portfolio.^Gamma;
        Exp_U=mean(U);
        CE(n,m)=Exp_U.^(1/Gamma);
              
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ALPHA_1s,ALPHA_2s] = meshgrid(Alpha_1s,Alpha_2s);
figure
hf=surf(ALPHA_1s,ALPHA_2s,CE');

% plot plane for radial allocation
a1=[Alpha_1s(1) Alpha_1s(1) Alpha_1s(end) Alpha_1s(end)];
a2=Line*a1;
Zz=[0 .9*max(max(CE)) .9*max(max(CE)) 0];
hold on
h=fill3(a1,a2,Zz,Zz);
grid on