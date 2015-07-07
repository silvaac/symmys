clc; clear; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameters

Mu=[0 0]';
s=[1 1];
rhos=[-.99 : .01 : .99];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n=1:length(rhos)
    rho=rhos(n);
    Sigma=[s(1)^2     rho*s(1)*s(2)
        rho*s(1)*s(2)    s(2)^2];

    [Expected_Value,Covariance,Standard_Deviation,Correlation]=LogNormalParam2Statistics(Mu,Sigma);
    Lambda=eig(Covariance);

    Cs(n)=Correlation(1,2);
    CRs(n)=min(Lambda)/max(Lambda);
end

figure

subplot(2,1,1)
plot(rhos,Cs)
ylim([-1 1])
grid on
xlabel('\rho')
ylabel('correlation')

subplot(2,1,2)
plot(rhos,CRs)
ylim([0 1])
grid on
xlabel('\rho')
ylabel('condition ratio')