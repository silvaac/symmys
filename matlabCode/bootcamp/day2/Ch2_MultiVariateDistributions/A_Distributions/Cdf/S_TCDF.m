% this script computes the cdf of a bivariate t distribution at a grid of points
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clear; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters
nu=20;
m=[0
    0];
s=[1
    1];
r=-.9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate grid according to natural disribution's limits
C=[1 r;
    r 1];
S=diag(s)*C*diag(s);

NumSimul=100;
Ones=ones(NumSimul,1);
Y=Ones*m' + (Ones*s').*mvtrnd(C,nu,NumSimul);
Grid1=[min(Y(:,1))  :  (max(Y(:,1))-min(Y(:,1)))/10  : max(Y(:,1))];
Grid2=[min(Y(:,2))  :  (max(Y(:,2))-min(Y(:,2)))/10  : max(Y(:,2))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute CDF on grid
IntegLim=abs(tinv(10^(-6),nu));
J=length(Grid1);
K=length(Grid2);
for j=1:J
    Countdown=J-j+1
    for k=1:K
        y=[Grid1(j)
            Grid2(k)];
        x=(y-m)./s;
        xmin = - IntegLim;
        xmax = min(x(1),IntegLim);
        ymin = - IntegLim;
        ymax = min(x(2),IntegLim);
        tol=10^(-6);
        method='quad';
        F(j,k) = dblquad('BivStandardTPDF',xmin,xmax,ymin,ymax,tol,method,r,nu);

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
[X,Y]=meshgrid(Grid1,Grid2);
surf(X,Y,F')