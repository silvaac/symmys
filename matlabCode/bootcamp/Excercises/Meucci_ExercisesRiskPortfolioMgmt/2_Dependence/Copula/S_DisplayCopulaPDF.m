clc; clear;
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters
Mu=[0  0]';     
r=.0;            
sigmas=[1 2]';    
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

         f_U(j,k)=TCopulaPDF(u,nu,Mu,Sigma);

    end
end
    
[G1,G2]=meshgrid(GridSide1,GridSide2);
figure
surf(G1,G2,f_U')
zlim([0 3])
xlabel('U_1')
ylabel('U_2')
zlabel('Copula pdf')
