% this script illustrates PCA dimension reduction in a tri-variate case
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; close all; clear;

Mu=[0.1 0.2 0.1];
Sig=[.25 .15 .2];
rho_12=.3;  rho_13=.3;  rho_23=.3;   

NumSimulations=10000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sigma=diag(Sig)*[1 rho_12 rho_13;  rho_12 1 rho_23;  rho_13 rho_23 1]*diag(Sig);

% compute sample: the market X is the first entry, and the two factors F1 and F2 are the last two entries
Y=mvnrnd(Mu,Sigma,NumSimulations);
Simul_X=Y;%exp(Y);

% compute recovered variables
Expected_Value=mean(Simul_X)';
Covariance=cov(Simul_X);
[EigenVectors,EigenValues] = pcacov(Covariance);
E_k=EigenVectors(:,1:2);
Recovered_X=[];
for t=1:NumSimulations
    X=Simul_X(t,:)';
    X_p=Expected_Value+E_k*E_k'*(X-Expected_Value);
   
    Recovered_X=[Recovered_X
        X_p' ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute ellipsoid
Scale=3;
Theta = [pi/2 : pi/50 : pi];
Phi = [0 : pi/80 : 2*pi];
Ellipsoid_x=zeros(length(Theta),length(Phi)); 
Ellipsoid_y=Ellipsoid_x; 
Ellipsoid_z=Ellipsoid_x;
for i=1:length(Theta)
    for j=1:length(Phi)
        
        y=[sin(Theta(i))*cos(Phi(j)); sin(Theta(i))*sin(Phi(j)); cos(Theta(i))];
        Ellipsoid = Expected_Value + Scale*EigenVectors*diag(sqrt(EigenValues))*y;
        
        Ellipsoid_x(i,j)= Ellipsoid(1);
        Ellipsoid_y(i,j)= Ellipsoid(2);
        Ellipsoid_z(i,j)= Ellipsoid(3);
        
    end
end

% compute plane of first two principal axes
Range=[-Scale*1.1 :Scale/10: Scale*1.1];
PrincipalPlane_X=zeros(length(Range),length(Range)); 
PrincipalPlane_Y=PrincipalPlane_X; 
PrincipalPlane_Z=PrincipalPlane_X;
for i=1:length(Range)
    for j=1:length(Range)
        VectorInPlane = Expected_Value + sqrt(EigenValues(1))*EigenVectors(:,1)*Range(i)...
                        + sqrt(EigenValues(2))*EigenVectors(:,2)*Range(j) ...
                    - sqrt(EigenValues(3))*EigenVectors(:,3)*.02; % cheat for cosmetics
        PrincipalPlane_X(i,j) = VectorInPlane(1);
        PrincipalPlane_Y(i,j) = VectorInPlane(2);
        PrincipalPlane_Z(i,j) = VectorInPlane(3);    
        
    end
end

% plot original sample
figure
h=plot3(Simul_X(:,1),Simul_X(:,2),Simul_X(:,3),'.');
% plot the principal plane
hold on
h=mesh(PrincipalPlane_X,PrincipalPlane_Y,PrincipalPlane_Z);
% plot the ellipsoid
hold on
h=contour3(Ellipsoid_x,Ellipsoid_y,Ellipsoid_z,200);
colormap('gray')
axis equal
grid on

% plot recovered variables
figure
h=plot3(Recovered_X(:,1),Recovered_X(:,2),Recovered_X(:,3),'.');
% plot the principal plane
hold on
h=mesh(PrincipalPlane_X,PrincipalPlane_Y,PrincipalPlane_Z);
% plot the ellipsoid
hold on
h=contour3(Ellipsoid_x,Ellipsoid_y,Ellipsoid_z,200);
colormap('gray')
axis equal
grid on