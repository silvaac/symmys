% this script illustrates regression dimension reduction in a simple bivariate case
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; close all; clear;

% the market X corresponds to the first entry, and the two factors F1 and F2 to the last two entries
Mu=[0.1 0.2 0.1];
Sig=[.25 .15 .25];
rho_F1X=.7;  rho_F2X=.7;  rho_F1F2=.5;   

NumSimulations=10000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sigma=diag(Sig)*[1 rho_F1X rho_F2X;  rho_F1X 1 rho_F1F2;  rho_F2X rho_F1F2 1]*diag(Sig);

% compute sample: the market X is the first entry, and the two factors F1 and F2 are the last two entries
Y=mvnrnd(Mu,Sigma,NumSimulations);
Simul_XFF=Y;%exp(Y);

% compute recovered variables
Expected_Value=mean(Simul_XFF)';
Covariance=cov(Simul_XFF);

ExpVal_X=Expected_Value(1);
Covariance_XX=Covariance(1,1);
ExpVal_F=Expected_Value(2:3);
Covariance_FF=Covariance(2:3,2:3);
Covariance_XF=Covariance(1,2:3);

B=Covariance_XF*inv(Covariance_FF);

Recovered_XFF=[];
for t=1:NumSimulations
    X=Simul_XFF(t,2);
    Y=Simul_XFF(t,3);
    Z=ExpVal_X+B*([X;Y]-ExpVal_F);
    Recovered_XFF=[Recovered_XFF
        Z X Y ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots

% compute regression plane 
Range=[-3 : .5 : 3];
RegressionPlane_X=zeros(length(Range),length(Range)); 
RegressionPlane_Y=RegressionPlane_X; 
RegressionPlane_Z=RegressionPlane_X;
Plane_X=zeros(length(Range),length(Range)); Plane_Y=RegressionPlane_X; Plane_Z=RegressionPlane_X;
MaxStdev=max( sqrt(Covariance_FF(1,1)), sqrt(Covariance_FF(2,2)));
for i=1:length(Range)
    for j=1:length(Range)
        x=ExpVal_F(1)+MaxStdev*Range(i);
        y=ExpVal_F(2)+MaxStdev*Range(j);
        z=ExpVal_X+B*([x;y]-ExpVal_F);
        
        RegressionPlane_X(i,j) = x;
        RegressionPlane_Y(i,j) = y;
        RegressionPlane_Z(i,j) = z;
    end
end


figure % original variables
% plot the regression plane
hl=surf(RegressionPlane_X,RegressionPlane_Y,RegressionPlane_Z);
% plot random simulations 
hold on
h=plot3(Simul_XFF(:,2),Simul_XFF(:,3),Simul_XFF(:,1),'.');
colormap('white')
grid on

figure % recovered variables
% plot the regression plane
hl=surf(RegressionPlane_X,RegressionPlane_Y,RegressionPlane_Z);
% plot the recovered variables
hold on
h=plot3(Recovered_XFF(:,2),Recovered_XFF(:,3),Recovered_XFF(:,1),'.');
colormap('white')
grid on