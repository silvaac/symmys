% this script shows that the eigenvectors of a Toplitz matrix have a Fourier basis structure
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

close all; clear; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=200; % dimension of the matrix
Decay=.9;      % decay factor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T=eye(N);
for n=1:N-1
    T=T+Decay^n*(diag(ones(N-n,1),n)+diag(ones(N-n,1),-n));
end
[EigVect,EigVals]=eig(T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for n=N-1:N
    hold on
    h=plot(EigVect(:,n));
    set(h,'color',[rand() rand() rand()]);
end
grid on