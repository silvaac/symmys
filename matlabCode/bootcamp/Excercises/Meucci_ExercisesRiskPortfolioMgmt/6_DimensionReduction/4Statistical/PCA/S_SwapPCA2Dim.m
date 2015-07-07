% this script performs the principal component analysis of a simplified two-point swap curve. 
% it computes and plots, among others, 
% 1. the invariants, namely rate changes
% 2. the location-dispersion ellipsoid of rates along with the 2-d location-dispersion ellipsoid
% 3. the effect on the curve of the two uncorrelated principal factors 
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci


clear; clc; close all;
load DB_Swap2y4y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% current curve
Current_Curve=Rates(end,:)
figure
h=plot([2 4],Current_Curve);
title('Current curve')
xlabel('Time to maturity, years')
ylabel('Par swap rate, %')

% determine weekly invariants (changes in rates)
Keep=[1:5:length(Dates)];

Rates=Rates(Keep,:);
X=Rates(2:end,:)-Rates(1:end-1,:);

Dates=Dates(Keep);
Dates(1)=[];

RunAnalysis(Dates,X(:,1),'weekly 2yr rates');
RunAnalysis(Dates,X(:,2),'weekly 4yr rates');

% scatter plot of  invariants
figure
plot(X(:,1),X(:,2),'.')
m=0*mean(X)'; % estimator shrunk to zero
S=cov(X);
TwoDimEllipsoid(m,S,2,1,0);
xlabel('2yr rate')
ylabel('4yr rate')

% perform PCA
[EVecs,EVals]=eig(S);
EVals=diag(EVals);
% sort eigenvalues in decreasing order
[dummy,Index]=sort(-EVals);
EigVals=EVals(Index);
EigVecs=EVecs(:,Index);

% plot eigenvectors
figure
h=plot([2 4],EigVecs(:,1),'r');
hold on 
h=plot([2 4],EigVecs(:,2),'g');
legend('1st factor loading','2nd factor loading');
xlabel('Time to maturity, years')

% factors
F=X*EigVecs;
F_std=std(F);

figure % 1-st factor effect
h=plot([2 4],Current_Curve);
hold on 
h=plot([2 4],Current_Curve'+F_std(1)*EigVecs(:,1),'r');
hold on 
h=plot([2 4],Current_Curve'-F_std(1)*EigVecs(:,1),'g');
legend('base','+1 sd of 1st fact','-1 sd of 1st fact')
xlabel('Time to maturity, years')

figure % 2-nd factor effect
h=plot([2 4],Current_Curve);
hold on 
h=plot([2 4],Current_Curve'+F_std(2)*EigVecs(:,2),'r');
hold on 
h=plot([2 4],Current_Curve'-F_std(2)*EigVecs(:,2),'g');
legend('base','+1 sd of 2nd fact','-1 sd of 2nd fact')
xlabel('Time to maturity, years')

% generalized R2
R2=cumsum(EigVals)/sum(EigVals)  % first entry: one factor, second entry: both factors