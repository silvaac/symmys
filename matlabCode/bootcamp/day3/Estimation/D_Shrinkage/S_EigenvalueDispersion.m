% this script displays the scattering-of-eigenvalues phenomenon 
% see Sec. 4.4 in "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clc; clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
N = 40;
Mu = 0*ones(N,1);
Sigma=eye(N);

NumSimulations=30;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sample=[round(1*N) : 50 : 20*N];

% compute true eigenvalues
[EVec,EVal]=eig(Sigma);
[dummy,Index]=sort(-diag(EVal));
EVec=EVec(:,Index);
EVal=diag(EVal(Index,Index));

% compute eigenvalues of sample estimator
Store_EVal_Hat=[];
for i=1:length(Sample)
    CyclesToGo=length(Sample)-i+1
    T=Sample(i);
    EVal_Hat=0;
    for n=1:NumSimulations
        X=mvnrnd(Mu,Sigma,T);
        Sigma_Hat=cov(X);

        [EVec_Hat,L]=eig(Sigma_Hat);
        [dummy,Index]=sort(-diag(L));
        L=diag(L(Index,Index));

        EVal_Hat=EVal_Hat+L;
    end
    EVal_Hat=EVal_Hat/NumSimulations;

    Store_EVal_Hat=[Store_EVal_Hat
        EVal_Hat'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots

figure
% plot eigenvalues of sample estimator
[X,Y]=meshgrid([1:N],Sample);
h1=surf(X,Y,Store_EVal_Hat);

% plot true eigenvalues
hold on
x=[1:N];
y=Sample(end)+0*EVal;
z=EVal;
h=plot3(x,y,z);

grid on
