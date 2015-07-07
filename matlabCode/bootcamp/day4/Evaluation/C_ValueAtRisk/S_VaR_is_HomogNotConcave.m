% this script shows that the quantile-based index of satisfaction is not necessarily concave
% i.e. the VaR is not necessarily convex
% however, the VaR is homogeneous and thus the decomposition of VaR in its contributors applies
% see "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clear;  close all;  clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs (discrete market)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% market parameters
Rho=-.99;
df=10;
G=20;

NumSimulations=10000;

Confidence=.3;

Line=2/3; % radial line for homogeneity check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sigma=[1 Rho; Rho 1];
A = chol(Sigma);
R_1 = unidrnd(G,NumSimulations/2,1)-(G+1)/2;
R_1 = [R_1
    -R_1];
R_2=mvtrnd(1,df,NumSimulations/2);
R_2 = [R_2
    - R_2];
R=[R_1 R_2]*A;
Prices=exp(R);

Alpha_1s=[0:.05:1];
Alpha_2s=[0:.05:1];
VaR=zeros(length(Alpha_1s),length(Alpha_2s));
Index=round(Confidence*NumSimulations);
for n=1:length(Alpha_1s)
    Countdown = length(Alpha_1s)-n
    for m=1:length(Alpha_2s)
        Alpha=[Alpha_1s(n); Alpha_2s(m)];
        Portfolio=Prices*Alpha;
        Sort_Portfolio=sort(Portfolio);
        VaR(n,m)=Sort_Portfolio(Index);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ALPHA_1s,ALPHA_2s] = meshgrid(Alpha_1s,Alpha_2s);

figure
hf=surf(ALPHA_1s,ALPHA_2s,VaR');

% plot plane for radial allocation
a1=[Alpha_1s(1) Alpha_1s(1) Alpha_1s(end) Alpha_1s(end)];
a2=Line*a1;
Zz=[0 .9*max(max(VaR)) .9*max(max(VaR)) 0];
hold on
h=fill3(a1,a2,Zz,Zz);
grid on