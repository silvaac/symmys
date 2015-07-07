% this script simulates a matrix-variate normal distribution
% it computes the summary statistics both in simulations and analytically
% see Sec. 2.6.2 in "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; clear; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters

% location
M=[2   -1      
   2    5
   4   -3];

% dispersion
S=[2   .5
   .5   2];
Sigma=[1     -.2     0
       -.2    4     -.7
       0     -.7    2];

% pick colums to display cross-covariances   
Col_j=2;
Col_k=1;

% pick rows to display cross-covariances   
Row_m=2;
Row_n=3;

NumSimulations=100000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rearrange data
N=size(M,1);
K=size(M,2);
Vec_M=reshape(M,N*K,1);
Kron_S=kron(S,Sigma);

% generate sample
Vec_X = mvnrnd(Vec_M,Kron_S,NumSimulations);
X=zeros(NumSimulations,N,K);
for n=1:N
    for k=1:K
        X(:,n,k)=Vec_X(:,(k-1)*N+n  );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute summary statistics
ExpVal=M
SampleMean=squeeze(mean(X,1))


SampleCov_ColjColk=zeros(N,N);
for m=1:N
    for n=1:N
        SS=cov(X(:,m,Col_j),X(:,n,Col_k));
        SampleCov_ColjColk(m,n)=SS(1,2);
    end
end
SampleCov_ColjColk
Cov_ColjColk=S(Col_j,Col_k)*Sigma

SampleCov_RowmRown=zeros(K,K);
for j=1:K
    for k=1:K
        SS=cov(X(:,Row_m,j),X(:,Row_n,k));
        SampleCov_RowmRown(j,k)=SS(1,2);
    end
end
SampleCov_RowmRown
Cov_RowmRown=Sigma(Row_m,Row_n)*S