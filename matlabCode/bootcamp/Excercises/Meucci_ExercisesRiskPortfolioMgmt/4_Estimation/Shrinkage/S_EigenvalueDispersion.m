% this script display the sample eigenvalues dispersion phenomenon
% see "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clc; clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
N = 50;
SampleLenght=[N : N : 10*N];

NumSimulations=50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mu = zeros(N,1);
Sigma= eye(N);

% compute true eigenvalues
[EVec,EVal]=eig(Sigma);
[dummy,Index]=sort(-diag(EVal));
EVec=EVec(:,Index);
EVal=diag(EVal(Index,Index));

% compute eigenvalues of sample estimator
Store_EVal_Hat=[];
for i=1:length(SampleLenght)
   
    T=SampleLenght(i);
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
[X,Y]=meshgrid([1:N],SampleLenght/N);
h1=surf(X,Y,Store_EVal_Hat);
xlabel('eigenvalue #')
ylabel('sample lenght/N')
grid on
