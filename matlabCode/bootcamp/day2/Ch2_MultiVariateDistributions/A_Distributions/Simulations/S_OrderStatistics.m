% this script shows that the pdf of the r-th order statistics is concentrated 
% around the quantile wiht confidence r/T
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci


close all; clc; clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input  
mu=0;
s=1;
nu=120;

T=100;
r=40;

NumSimul=30000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=mu+s*trnd(nu,NumSimul,T);

X_T=sort(X,2);
X_rT=X_T(:,r);

q=mu+s*tinv(r/T,nu);

F=tcdf((X_rT-mu)/s,nu);
f=1/s*tpdf((X_rT-mu)/s,nu);
pdf_rT = gamma(T+1)/(gamma(r)*gamma(T-r+1))*(F.^(r-1)).*((1-F).^(T-r)).*f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
NumBins=round(10*log(NumSimul));
hist(X_rT,NumBins)
[n,D]=hist(X_rT,NumBins);
RescaleFactor=(D(2)-D(1))*NumSimul;
RescaledPdf_rT=pdf_rT *RescaleFactor;
hold on
h=plot(X_rT,RescaledPdf_rT,'.');
set(h,'color','r')
hold on
h=plot(q,0,'.');
set(h,'color','r','markersize',20)
grid on