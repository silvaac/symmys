% this script simulates the copula of the diagonal elements of a 2x2 Wishart distribution
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clear;  close all;  clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=[1 2];
r=0.6;
nu=20;
NumSimul=10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sigma=diag(s)*[1 r;r 1]*diag(s);


W_xx=[]; W_yy=[]; 
for j=1:NumSimul
  X = mvnrnd(zeros(2,1),Sigma,nu);
  W = X'*X;
   
  W_xx=[W_xx
    W(1,1)];
  W_yy=[W_yy
    W(2,2)];
 
end

a_xx=nu/2;
b_xx=2*Sigma(1,1);
U_xx=gamcdf(W_xx,a_xx,b_xx); % grade 1 simulation
a_yy=nu/2;
b_yy=2*Sigma(2,2);
U_yy=gamcdf(W_yy,a_yy,b_yy); % grade 2 simulation
Copula=[U_xx U_yy];			      % copula

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