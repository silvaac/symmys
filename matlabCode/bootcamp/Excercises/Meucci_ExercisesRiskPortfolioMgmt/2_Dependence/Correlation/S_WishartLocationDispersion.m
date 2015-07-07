% this script computes the location-dispersion ellipsoid of the 
% normalized (unit variance, zero expectation) first diagonal and
% off-diagonal elements of a 2x2 Wishart distribution
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clear;  close all;  clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=[1 1];
r=-0.9;
nu=5;
Num_Simluations=10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sigma=diag(s)*[1 r;r 1]*diag(s);

W_xx=[]; W_yy=[]; W_xy=[];
Vec_W=[];
Dets=[];
Traces=[];
for j=1:Num_Simluations
  X = mvnrnd(zeros(2,1),Sigma,nu);
  W = X'*X;
 
  Dets=[Dets det(W)];
  Traces=[Traces trace(W)];
  
  W_xx=[W_xx
    W(1,1)];
  W_yy=[W_yy
    W(2,2)];
  W_xy=[W_xy
    W(1,2)];

  Vec_W=[Vec_W
      reshape(W,1,4)];
end

% compute expected values of W_xx and W_xy, see (2.227) in "Risk and Asset Allocation - Springer
E_xx=nu*Sigma(1,1);
E_xy=nu*Sigma(1,2);

% compute covariance matrix of W_xx and W_xy, see (2.228) in "Risk and Asset Allocation - Springer
m=1;n=1;p=1;q=1;
var_Wxx=nu*(Sigma(m,p)*Sigma(n,q)+Sigma(m,q)*Sigma(n,p));
m=1;n=2;p=1;q=2;
var_Wxy=nu*(Sigma(m,p)*Sigma(n,q)+Sigma(m,q)*Sigma(n,p));
m=1;n=1;p=1;q=2;
cov_Wxx_Wxy=nu*(Sigma(m,p)*Sigma(n,q)+Sigma(m,q)*Sigma(n,p));

S_xx_xy=[var_Wxx cov_Wxx_Wxy 
    cov_Wxx_Wxy var_Wxy];

% compute X_1 and X_2, i.e. normalized version of W_xx and W_xy
X_1=(W_xx-E_xx)/sqrt(var_Wxx);
X_2=(W_xy-E_xy)/sqrt(var_Wxy);
X=[X_1 X_2];

% compute expected value and covariance of X_1 and X_2
E=[0 
    0];
E_hat=mean(X)'

S=diag(1./[sqrt(var_Wxx) sqrt(var_Wxy)]) *  S_xx_xy   * diag(1./[sqrt(var_Wxx) sqrt(var_Wxy)])
S_hat=cov(X)

figure
plot(X_1,X_2,'.')
Scale=1;
PlotEigVectors=1;
PlotSquare=0;
TwoDimEllipsoid(E,S,Scale,PlotEigVectors,PlotSquare)
xlabel('X_1')
ylabel('X_2')