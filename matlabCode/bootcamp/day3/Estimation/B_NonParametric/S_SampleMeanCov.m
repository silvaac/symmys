% this script shows that the sample mean and sample covariance correspond
% to the ellipsoid that best fits the data (i.e. least volume) among all those who fit the
% data (i.e. average Mahlanobis distance = 1)
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; close all; clear;

Mu=[1 1]';
Sigma=[1 1.5];
Rho=-.5;

T=100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%process inputs
Correlations=[1 Rho; Rho 1];
Covariance=diag(Sigma)*Correlations*diag(Sigma);

% compute sample mean and covariance
X=mvnrnd(Mu,Covariance,T);
Sample_Mean=mean(X)';
Sample_Covariance=cov(X);

% compute random location and dispersion parameters such that the average square Mahalanobis distance be one
m=Sample_Mean+randn(2,1);
A=rand(2,2)-.5; 
S=A*A'; 
Total_Mahalanobis=0;
inv_S=inv(S);
for t=1:T 
    x=X(t,:)';
    Total_Mahalanobis = Total_Mahalanobis + 1/T*(x-m)'*inv_S*(x-m);
end
S=S*Total_Mahalanobis;  % constraint: 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
figure
h=plot(X(:,1),X(:,2),'.');
hold on
h1=TwoDimEllipsoid(Sample_Mean,2*Sample_Covariance,1,0,0);
set(h1,'color','r','linewidth',2)

hold on
h2=TwoDimEllipsoid(m,S,1,0,0);
legend([h1 h2],'sample mean/cov','random')