% this script displays the pdfs of the log-t distribution
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci
% formula (2.213)

clc; clear; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters
Mu=[1  0]';     
r=.7;            
sigmas=[1.3 1.3]';    
nu=3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=[1 r;r 1];
Sigma=diag(sigmas)*C*diag(sigmas);

% generate sample to define grid
NumSimul=10000;
Z=mvtrnd(C,nu,NumSimul);
Y=ones(NumSimul,1)*Mu'+Z*diag(sigmas);  
X=exp(Y);

Percentile=.2;
Max=prctile(X(:,1),100*(1-Percentile));
Min=prctile(X(:,1),100*Percentile);
Step=(Max-Min)/50;
GridSide1=[Min : Step : Max];
Max=prctile(X(:,2),100*(1-Percentile));
Min=prctile(X(:,2),100*Percentile);
Step=(Max-Min)/50;
GridSide2=[Min : Step : Max];

% compute pdf on grid
for j=1:length(GridSide1)
    for k=1:length(GridSide2)
        x=[GridSide1(j)
            GridSide2(k)];
        f(j,k)=LogTPDF(x,nu,Mu,Sigma);
    end
end

% display pdf
[G1,G2]=meshgrid(GridSide1,GridSide2);
figure
surf(G1,G2,f')