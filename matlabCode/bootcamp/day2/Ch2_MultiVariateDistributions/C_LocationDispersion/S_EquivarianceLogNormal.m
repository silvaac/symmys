% this script shows the affine equivariance properties of expected value
% and covariance in a lognormal example
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; clear; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters
NumSimulations=50000;

% distribution
Mu=[0.1 0.1]';  
s=[0.3 0.2];
r=-0.5;

% invertible affine transformation
m=[.5 .4]';
B=[-1 2;.5 -2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate lognormal sample
Sigma=[s(1)^2     r*s(1)*s(2)
    r*s(1)*s(2)    s(2)^2];

Z = mvnrnd(Mu,Sigma,NumSimulations);
X = exp(Z);
Y = ones(NumSimulations,1)*m'+X*B';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute location and dispersion (analytical)
Expected_Value_X=exp(Mu+(1/2)*diag(Sigma));
Covariance_X=exp(Mu+(1/2)*diag(Sigma))*exp(Mu+(1/2)*diag(Sigma))'.*(exp(Sigma)-1);
Expected_Value_Y=m+B*Expected_Value_X;
Covariance_Y=B*Covariance_X*B';

% compute Mahalanobis distance
u=X-ones(NumSimulations,1)*Expected_Value_X';
ZX=sqrt(   sum( (u*inv(Covariance_X)).*u,2)   );

v=Y-ones(NumSimulations,1)*Expected_Value_Y';
ZY=sqrt(   sum( (v*inv(Covariance_Y)).*v,2)   );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
figure 

Scale=2;
PlotEigVectors=1;
PlotSquare=1;

h=plot(X(:,1),X(:,2),'.');
hold on 
TwoDimEllipsoid(Expected_Value_X,Covariance_X,Scale,PlotEigVectors,PlotSquare)
hold on 
h=plot(Y(:,1),Y(:,2),'.');
set(h,'color','g')
hold on 
TwoDimEllipsoid(Expected_Value_Y,Covariance_Y,Scale,PlotEigVectors,PlotSquare)
grid on

figure
NumBins=round(10*log(NumSimulations));

subplot(2,1,1)
hist(ZX,NumBins)
xlabel('Mahalanobis distance of X')
grid on
subplot(2,1,2)
hist(ZY,NumBins)
xlabel('Mahalanobis distance of Y=a+BX')
grid on