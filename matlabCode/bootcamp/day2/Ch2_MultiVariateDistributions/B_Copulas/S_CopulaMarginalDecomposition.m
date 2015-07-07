% this script generates a joint distribution with normal copula and arbitrary marginals
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

close all; clc; clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters
NumSimul=10000;

% input for bivariate normal distribution
NormCorr=-.8;
NormStDev=[1 3]';      % NOTE: this input plays no role in the final output
NormExpVal=[-2 5]';    % NOTE: this input plays no role in the final output

% input for first marginal
nu_1=9;
mu_1=-1;
sigmasq_1=2;

% input for second marginal
nu_2=7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate bivariate normal distribution
NormCorrMatrix = [1 NormCorr;NormCorr 1];
NormCovMatrix = diag(NormStDev)*NormCorrMatrix*diag(NormStDev);
Z = mvnrnd(NormExpVal,NormCovMatrix,NumSimul);

Z_1=Z(:,1);
Z_2=Z(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate copula
U_1 = normcdf(Z(:,1),NormExpVal(1),NormStDev(1));  % grade 1
U_2 = normcdf(Z(:,2),NormExpVal(2),NormStDev(2));  % grade 2
U=[U_1 U_2];   % joint realizations from the required copula

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate joint distribution
X_1 = mu_1+sqrt(sigmasq_1)*tinv(U_1,nu_1);
X_2 = chi2inv(U_2,nu_2);
X=[X_1 X_2];   % joint realizations from the required distribution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot joint distribution
figure % marginals: as expected, the histograms (pdf's) do NOT change as NormCorr varies
NumBins=round(10*log(NumSimul));
subplot(2,1,1)  % t distribution
hist(X_1,NumBins)
xlabel('t distribution');
grid on
subplot(2,1,2) % chi-square distribution
hist(X_2,NumBins)
xlabel('chi square distribution');
grid on

figure % joint
plot(X_1,X_2,'.')
grid on
xlabel('t distribution');
ylabel('chi square distribution');

break
figure % pdf
NumBins2D=round(sqrt(NumSimul)/5);
NumBins2D=[NumBins2D NumBins2D];
hist3(X,NumBins2D);
title('pdf joint distribution')
xlabel('t distribution');
ylabel('chi square distribution');
