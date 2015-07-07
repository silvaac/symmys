% this script shows that the pdf of the r-th order statistics is concentrated 
% around the quantile wiht confidence r/T
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

close all; clc; clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input  
mu=0.2;
s=0.25;
T=70;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pdf of r-th order statistic concentrated around the r/T quantile

rs=[1: T];
x=[0 : .01 : 2.5*exp(mu+s*s/2)];
F=logncdf(x,mu,s);
f=lognpdf(x,mu,s);
for n=1:length(rs)
    r=rs(n);
    
    pdf_rT = gamma(T+1)/(gamma(r)*gamma(T-r+1))*(F.^(r-1)).*((1-F).^(T-r)).*f;
    q=logninv(r/T,mu,s);

    hold on
    plot3(x,r/T+0*x,pdf_rT);
    hold on
    plot3(q,r/T,0,'.')
end
view([-20,60])
xlabel('x')
ylabel('r/T')
zlabel('pdf')
