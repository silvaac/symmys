% this script shows that the pdf of the r-th order statistics is concentrated 
% around the quantile wiht confidence r/T
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci
% formula (2.248) and (2.253)

close all; clc; clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input  
mu=0;
s=1;
nu=10;
T=70;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pdf of r-th order statistic concentrated around the r/T quantile

rs=[1: T];
x=mu + s*[-4: .01 : 4];
F=tcdf((x-mu)/s,nu);
f=1/s*tpdf((x-mu)/s,nu);
for n=1:length(rs)
    r=rs(n);
    
    pdf_rT = gamma(T+1)/(gamma(r)*gamma(T-r+1))*(F.^(r-1)).*((1-F).^(T-r)).*f;
    q=mu+s*tinv(r/T,nu);

    hold on
    plot3(x,r/T+0*x,pdf_rT);
    hold on
    plot3(q,r/T,0,'.')
end

