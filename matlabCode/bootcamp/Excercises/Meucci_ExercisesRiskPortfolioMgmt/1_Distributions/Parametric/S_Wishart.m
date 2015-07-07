% this script generates a sample from the 2x2 Wishart distribution 
% it shows that determinant and trace are positive, i.e. the matrix is positive
% it shows that the marginal diagonal are gamma-distributed
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clear;  close all;  clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=[1 1];
r=0.3;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute summary statistics (analytical and sample-based)
Expected_Value=reshape(nu*Sigma,1,4)
K_22=[1	 0	0	0
    0	0	1	0
    0	1	0	0
    0	0	0	1];
Covariance=nu*(eye(4)+K_22)*kron(Sigma,Sigma)

Sample_Mean=mean(Vec_W)
Sample_Covariance=cov(Vec_W)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots

figure % tri-variate joint
h=plot3(W_xx,W_xy,W_yy,'.');
xlabel('xx')
ylabel('xy')
zlabel('yy')
grid on

figure % bi-variate marginals
subplot(2,2,1)
plot(W_xx,W_xy,'.')
xlabel('xx')
ylabel('xy')
grid on
subplot(2,2,2)
plot(W_yy,W_xy,'.')
xlabel('yy')
ylabel('xy')
grid on
subplot(2,2,3)
plot(W_xx,W_yy,'.')
xlabel('xx')
ylabel('yy')
grid on

figure % marginals

NumBins=round(10*log(Num_Simluations));
subplot(2,2,1)
hist(W_xx,NumBins);
[n,D]=hist(W_xx,NumBins);
RescalePdf=(D(2)-D(1))*Num_Simluations;
y=gampdf(W_xx,nu/2,2*Sigma(1,1))*RescalePdf;
hold on
h=plot(W_xx,y,'.');
grid on
xlabel('xx')

subplot(2,2,2)
hist(W_yy,NumBins)
[n,D]=hist(W_yy,NumBins);
RescalePdf=(D(2)-D(1))*Num_Simluations;
y=gampdf(W_yy,nu/2,2*Sigma(2,2))*RescalePdf;
hold on
h=plot(W_yy,y,'.');
grid on
xlabel('yy')

subplot(2,2,3)
hist(W_xy,NumBins)
grid on
xlabel('xy')

figure % check positivity
subplot(2,1,1)
hist(Traces,NumBins)
xlabel('trace')
grid on
subplot(2,1,2)
hist(Dets,NumBins)
xlabel('determinant')
grid on