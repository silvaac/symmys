% this script displays the cdf of the normal and t copulas 
% (same as lognormal and log-t copulas respectively)
% see Sec. 2.2.2 in "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; clear; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters
Mu=[1  -10]';     
r=.5;            
sigmas=[1 1]';    
nu=100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sigma=diag(sigmas)*[1 r;r 1]*diag(sigmas);

% compute and display results
GridSide1=[.05:.05:.95];
GridSide2=GridSide1;

for j=1:length(GridSide1)
    for k=1:length(GridSide2)
        
        u=[GridSide1(j) GridSide2(k)]';

         F_U(j,k)=NormalCopulaCDF(u,Mu,Sigma);    
         %F_U(j,k)=TCopulaCDF(u,nu,Mu,Sigma);  
    end
end
    
[G1,G2]=meshgrid(GridSide1,GridSide2);
figure
surf(G1,G2,F_U')