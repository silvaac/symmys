% this script compare PCA (implicit factors) with regression (explicit factors) 
% dimension reduction in a simple bi-variate case
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; close all; clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ExpectedValue = [1 1]';
Standard_Deviations=[1.7 1.5];
Rho=.7;
NumSimulations=200;
Scale=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute PCA analysis
Covariance=diag(Standard_Deviations)*[1 Rho; Rho 1]*diag(Standard_Deviations);
[EigenVectors,EigenValues] = pcacov(Covariance);
EigenVectors(:,1)=-EigenVectors(:,1);
Rsquare_PCA=EigenValues(1)/sum(EigenValues);

% compute the ellipsoid in the r plane, solution to  ((R-Mean)' * Covariance^-1 * (R-Mean) ) = Scale^2                                   
Centered_Ellipse=[]; 
Angle = [0 : pi/500 : 2*pi];
Steps = length(Angle);
for i=1:Steps 
  y=[cos(Angle(i))                                   % normalized variables (parametric representation of the ellipsoid)
    sin(Angle(i))];
  Centered_Ellipse=[Centered_Ellipse Scale*EigenVectors*diag(sqrt(EigenValues))*y];  
end
Ellipse= ExpectedValue*ones(1,Steps) + Centered_Ellipse;

% compute pricipal planes (one axis in this case)
t=[-1.2 : 0.01 : 1.2];
MainVariationTrasl_Line = ExpectedValue*(1+0*t) + Scale*sqrt(EigenValues(1))*EigenVectors(:,1)*t;

% compute regression hyperplanes (lines in this case)
Regression_x = [ExpectedValue(1) - Scale*1.5*sqrt(Covariance(1,1)) : 0.01 : ExpectedValue(1) + Scale*1.5*sqrt(Covariance(1,1))];
Regression_y = ExpectedValue(2) + Covariance(2,1)/Covariance(1,1) * (Regression_x-ExpectedValue(1));

% compute sample
Simul_X=mvnrnd(ExpectedValue,Covariance,NumSimulations);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots 

% regression analysis

figure
h=plot(Simul_X(:,1),Simul_X(:,2),'.'); % plot the simulations
set(h,'color','k')
hold on
h=plot(Ellipse(1,:),Ellipse(2,:)); % plot the ellipse
set(h,'color','g','linewidth',2)
hold on
h=plot(MainVariationTrasl_Line(1,:),MainVariationTrasl_Line(2,:)); % plot PCA line
set(h,'color','r','linewidth',2)
hold on 
h=plot(Regression_x,Regression_y); % plot regression line
set(h,'color','b','linewidth',2)
axis equal;
grid on
