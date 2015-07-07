% this script displays the pdf of the copula 
% of the normal, lognormal, t, and log-t distribution
% see Sec. 2.2.2 in "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; clear; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters
Mu=[0.04  0.05]';     
r=.0;            
sigmas=[.25 .3]';    

nu=1;

Sigma=diag(sigmas)*[1 r;r 1]*diag(sigmas);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute and display results

GridSide1=[.05:.05:.95];
GridSide2=GridSide1;

for j=1:length(GridSide1)
    for k=1:length(GridSide2)
        u=[GridSide1(j)
            GridSide2(k)];

         % !!!!!   choose below which pdf to display by commenting those that should
         % not be computed !!!!
            
         %f_U(j,k)=NormalCopulaPDF(u,Mu,Sigma);    
         %f_U(j,k)=LognormalCopulaPDF(u,Mu,Sigma);
         f_U(j,k)=TCopulaPDF(u,nu,Mu,Sigma);
         %f_U(j,k)=LogTCopulaPDF(u,nu,Mu,Sigma);
    end
end
    
[G1,G2]=meshgrid(GridSide1,GridSide2);
figure
surf(G1,G2,f_U')
zlim([0 3])
