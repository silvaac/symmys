% this script 
% see "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clear; clc; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input data
load ('db_SwapRatesUS');
ChooseRates=[1 3]; % 1=2yr; 2=5yr; 3=10yr

Y=[Data(:,ChooseRates(1)+1) Data(:,ChooseRates(2)+1)];
X=Y(2:end,:)-Y(1:end-1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% non-parametric

M_NoPar=mean(X)'
S_NoPar=cov(X)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MLE

[Nu,Mu,Sigma]=StudentMLE(X);
M_MLE=Mu
S_MLE=Sigma*Nu/(Nu-2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shrinkage

[M_Shr,S_Shr]=Shrinkage(X)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Robust
[M_Rob,S_Rob]=HubertM(X)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot results
figure 

h=plot(X(:,1),X(:,2),'.');
set(h,'color','k','markersize',3)

hold on
E1=TwoDimEllipsoid(M_NoPar,S_NoPar,2,0,0);
set(E1,'color','b')

hold on
E2=TwoDimEllipsoid(M_MLE,S_MLE,2,0,0);
set(E2,'color','r')

hold on
E3=TwoDimEllipsoid(M_Shr,S_Shr,2,0,0);
set(E3,'color','k')

hold on
E4=TwoDimEllipsoid(M_Rob,S_Rob,2,0,0);
set(E4,'color','y')

legend([E1 E2 E3 E4],'non parametric','MLE','shrinkage','robust','location','best')
title('daily swap rate changes')
xlabel(Fields(ChooseRates(1)+1).Name)
ylabel(Fields(ChooseRates(2)+1).Name)