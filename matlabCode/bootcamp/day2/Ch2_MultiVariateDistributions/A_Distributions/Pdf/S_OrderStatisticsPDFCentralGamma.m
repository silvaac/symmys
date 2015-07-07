% this script shows that the pdf of the r-th order statistics is concentrated 
% around the quantile wiht confidence r/T
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci
% formula (2.248) and (2.253)

close all; clc; clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input  
nu=10;
s2=1.2;
T=70;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pdf of r-th order statistic concentrated around the r/T quantile
a=nu/2;
b=2*s2;

rs=[1: T];
x=[0 : .01 : nu*s2+3*sqrt(2*nu)*s2];
F=gamcdf(x,a,b);
f=gampdf(x,a,b);
for n=1:length(rs)
    r=rs(n);
    
    pdf_rT = gamma(T+1)/(gamma(r)*gamma(T-r+1))*(F.^(r-1)).*((1-F).^(T-r)).*f;
    q=gaminv(r/T,a,b);

    hold on
    plot3(x,r/T+0*x,pdf_rT);
    hold on
    plot3(q,r/T,0,'.')
end

