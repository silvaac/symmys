% this script displays the pdfs of the 1x1 Wishart distribution, 
% i.e. the gamma distribution
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci
% formula (2.230)

clc; clear; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters
nu=10;            
sigma_square=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate sample to define grid
NumSimul=10000;
a=nu/2;
b=2*sigma_square;

% generate sample using matlab functions
W=gamrnd(a,b,NumSimul,1);

Percentile=.05;
Max=prctile(W,100*(1-Percentile));
Min=prctile(W,100*Percentile);
Step=(Max-Min)/50;
Grid=[Min : Step : Max];

% compute pdf on grid (and double check with gamma distribution)
for j=1:length(Grid)
        w=Grid(j);
        f_Wishart(j)=WishartPDF(w,nu,sigma_square);
        f_Gamma(j)=gampdf(w,a,b);
end

% display pdf
figure
plot(Grid,f_Wishart,'.')
hold on
plot(Grid,f_Gamma,'r')
legend('Wishart','gamma')
grid on