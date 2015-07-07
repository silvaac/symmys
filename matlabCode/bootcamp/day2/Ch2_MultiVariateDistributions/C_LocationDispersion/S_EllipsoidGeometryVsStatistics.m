% this script shows the relationship between geometry, algebra and statistics 
% in the location-dispersion ellipsoid using a lognormal example
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci


clc; clear; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters

Mu=[.05 .1]';
s=[.1 .1]';
nu=40;
r=-.9;

NumSimulations=10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate sample
C=[1 r; r 1];
Ones=ones(NumSimulations,1);
Y = Ones*Mu' + (Ones*s').*mvtrnd(C,nu,NumSimulations);
X = exp(Y);
m=mean(X)';
S=cov(X);

Theta=[0 : pi/100 : 2*pi];
for n=1:length(Theta)
    th=Theta(n);
    
    e=[cos(th)  % versor
        sin(th)];

    % evaluate standard deviation on a one-dim projection (versor)
    Z=X*e;      
    SDev(n)=std(Z);  
    
    % radius of ellipsoid
    Radius(n)=(e'*inv(S)*e)^(-1/2);   
end

% compute min and max standard deviation and respective versor
[dd,min_n]=min(SDev);
e_min= [cos(Theta(min_n))
        sin(Theta(min_n))];
s_min=SDev(min_n);

[dd,max_n]=max(SDev);
e_max=[cos(Theta(max_n))
        sin(Theta(max_n))];
s_max=SDev(max_n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure 
% scatter plot simulations
h=plot(X(:,1),X(:,2),'.');
hold on 

% plot ellipsoid
Scale=2;
PlotEigVectors=1;
PlotSquare=1;
TwoDimEllipsoid(m,S,Scale,PlotEigVectors,PlotSquare)
grid on

% plot special directions defined by the max-min versors
Center=mean(X)'*.7; % de-center plot of special directions for better display
v=Scale*[-1:.1:1]';
Ones=1+0*v;

v_min=Ones*Center'+v*s_min*e_min';
hold on
h=plot(v_min(:,1),v_min(:,2),'r');
v_max=Ones*Center'+v*s_max*e_max';
hold on 
h=plot(v_max(:,1),v_max(:,2),'r');
axis equal
grid on

% plot statistics versus geometry
figure
Scaled_Theta=Theta/pi;
h=plot(Scaled_Theta,SDev);  % plot standard deviation as function of direction
hold on 
h=plot(Scaled_Theta,Radius,'r'); % plot radius of ellipsoid as function of direction
xlim([Scaled_Theta(1) Scaled_Theta(end)])
grid on
xlabel('theta/pi')
legend('st.dev on projection','radius of ellipsoid')
