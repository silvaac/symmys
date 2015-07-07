% this script demonstrates the recursive ML estimation of the location and scatter
% parameters of a multivariate Student t distribution
% see "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clear; clc; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
load ('DBUsSwapRates');
ChooseRates=[1 2]; % 1=2yr; 2=5yr; 3=10yr

Y=[Series(1).Data Series(3).Data];
X=Y(2:end,:)-Y(1:end-1,:);

Nus=[3 100];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computations
Tolerance=10^(-10);
for q=1:length(Nus)
    [Estimate(q).Mu_hat,Estimate(q).Sigma_hat] = MleRecursionForT(X,Nus(q),Tolerance);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figures
figure
h=plot(X(:,1),X(:,2),'.');
set(h,'color','b','markersize',4)
for q=1:length(Nus)
    hold on
    M=Estimate(q).Mu_hat;
    S=Estimate(q).Sigma_hat*Nus(q)/(Nus(q)-2);
    dd=TwoDimEllipsoid(M,S,2,0,0);
    set(dd,'color',.7*[rand() rand() rand()])
end
xlim([-.4 .4])
ylim([-.4 .4])
xlabel(Series(ChooseRates(1)).Name)
ylabel(Series(ChooseRates(2)).Name)