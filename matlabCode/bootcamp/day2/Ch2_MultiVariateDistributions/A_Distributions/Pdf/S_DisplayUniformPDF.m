% this script displays the pdfs of the uniform distribution on an ellipsoid
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci
% formula (2.145)

clc; clear; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters
Mu=[1  0]';     
r=-.7;            
sigmas=[1.3 1.3]';    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sigma=diag(sigmas)*[1 r;r 1]*diag(sigmas);

% generate sample to define grid
NumSimul=10000;
X=mvnrnd(Mu,Sigma,NumSimul);

Percentile=.1;
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
        f(j,k)=UniformEllipticalPDF(x,Mu,Sigma);
    end
end

% display pdf
[G1,G2]=meshgrid(GridSide1,GridSide2);
figure
surf(G1,G2,f')