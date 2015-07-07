% this script simulates the copula of a bi-variate log-t distribution
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci


clc; clear; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters
Mu=[1  -1]';     
r=-.3;            
sigma=[2 10]';    
nu=200;

NumSimul=20000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=[1 r;r 1];
Std_X=mvtrnd(C,nu,NumSimul);  % standardized bi-variate t simulation
Y_1=Mu(1)+sigma(1)*Std_X(:,1);
Y_2=Mu(2)+sigma(2)*Std_X(:,2);
Y=[Y_1 Y_2]; % bi-variate t simulation

X=exp(Y);    % bi-variate log-t simulation

U_1=tcdf( (log(X(:,1))-Mu(1))/sigma(1) ,nu); % grade 1 simulation
U_2=tcdf( (log(X(:,2))-Mu(2))/sigma(2) ,nu); % grade 1 simulation
Copula=[U_1 U_2];                                 % copula

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure 
% marginals
NumBins=round(10*log(NumSimul));

subplot('Position',[.05 .3 .2 .6]) 
[n,D]=hist(Copula(:,2),NumBins);
barh(D,n,1);
[y_lim]=get(gca,'ylim')
set(gca,'xtick',[])
grid on

subplot('Position',[.3 .05 .6 .2]) 
[n,D]=hist(Copula(:,1),NumBins);
bar(D,n,1);
[x_lim]=get(gca,'xlim')
set(gca,'ytick',[])
grid on

% scatter plot
subplot('Position',[.3 .3 .6 .6]) 
h=plot(Copula(:,1),Copula(:,2),'.');
set(gca,'xlim',x_lim,'ylim',y_lim)
grid on
xlabel('grade 1');
ylabel('grade 2');

% 3-d histogram (~rescaled pdf)
NumBins3d=round(sqrt(NumSimul)/5);
figure
hist3(Copula(:,[1 2]),[NumBins3d NumBins3d]);