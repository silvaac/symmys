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
sigmasq_1=2;

mu_2=0;
sigmasq_2=.04;

% input for second marginal
nu_2=7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate bivariate normal distribution
NormCorrMatrix = [1 NormCorr;NormCorr 1];
NormCovMatrix = diag(NormStDev)*NormCorrMatrix*diag(NormStDev);
Z = mvnrnd(NormExpVal,NormCovMatrix,NumSimul);

Z_1=Z(:,1);
Z_2=Z(:,2);
% plots
figure % marginals: as expected, they are normal
NumBins=round(10*log(NumSimul));
subplot(2,1,1)  
hist(Z_1,NumBins)
xlabel('normal 1');
grid on
subplot(2,1,2) 
hist(Z_2,NumBins)
xlabel('normal 2');
grid on
figure % joint
plot(Z_1,Z_2,'.')
grid on
xlabel('normal 1');
ylabel('normal 2');

figure % pdf
NumBins2D=round(sqrt(100*log(NumSimul)));
NumBins2D=[NumBins2D NumBins2D];
hist3(Z,NumBins2D);
xlabel('normal 1')
xlabel('normal 2')
title('pdf normal')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate copula
U_1 = normcdf(Z(:,1),NormExpVal(1),NormStDev(1));  % grade 1
U_2 = normcdf(Z(:,2),NormExpVal(2),NormStDev(2));  % grade 2
U=[U_1 U_2];   % joint realizations from the required copula

% plot copula
figure % marginals: as expected, they are uniform
NumBins=round(10*log(NumSimul));
subplot(2,1,1)  
hist(U_1,NumBins)
xlabel('grade 1');
grid on
subplot(2,1,2) 
hist(U_2,NumBins)
xlabel('grade 2');
grid on

figure % joint
plot(U_1,U_2,'.')
grid on
xlabel('grade 1');
ylabel('grade 2');

figure % pdf
NumBins2D=round(sqrt(100*log(NumSimul)));
NumBins2D=[NumBins2D NumBins2D];
hist3(U,NumBins2D);
xlabel('grade 1');
ylabel('grade 2');
title('pdf copula')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate joint distribution

a=nu_1/2;
b=2*sigmasq_1;
X_1 = gaminv(U_1,a,b);

sigma_2=sqrt(sigmasq_2);
X_2 = logninv(U_2,mu_2,sigma_2);

X=[X_1 X_2];   % joint realizations from the required distribution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot joint distribution
figure % marginals: as expected, the histograms (pdf's) do NOT change as NormCorr varies
NumBins=round(10*log(NumSimul));
subplot(2,1,1)  % t distribution
hist(X_1,NumBins)
xlabel('gamma');
grid on
subplot(2,1,2) % chi-square distribution
hist(X_2,NumBins)
xlabel('lognormal');
grid on

figure % joint
plot(X_1,X_2,'.')
grid on
xlabel('gamma');
ylabel('lognormal');

figure % pdf
NumBins2D=round(sqrt(100*log(NumSimul)));
NumBins2D=[NumBins2D NumBins2D];
hist3(X,NumBins2D);
xlabel('gamma');
ylabel('lognormal');
title('pdf joint distribution')