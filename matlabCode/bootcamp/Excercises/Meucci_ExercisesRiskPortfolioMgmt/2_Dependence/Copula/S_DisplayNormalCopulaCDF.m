close all;
clc; clear;

Mu=[0  0]';       % exp values
r=-.99999;            % correlation
sigmas=[1 1]';     % st. deviations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sigma=diag(sigmas)*[1 r;r 1]*diag(sigmas);

GridSide1=[.05:.05:.95];
GridSide2=GridSide1;
NumGrid=length(GridSide1);

F_U=zeros(NumGrid);
for j=1:NumGrid
    for k=1:NumGrid;

        u=[GridSide1(j)
            GridSide2(k)];

        x=norminv(u,Mu,sigmas);

        F_U(j,k) = mvncdf(x,Mu,Sigma);

    end
end

[G1,G2]=meshgrid(GridSide1,GridSide2);
surf(G1,G2,F_U')
xlabel('U_1')
ylabel('U_2')
zlabel('Copula cdf')

