function PlotNIWMarginals(Mu_Simul,InvSigma_Simul,Mu_0,T_0,Sigma_0,Nu_0,Legend)
% this function plots numerically and analytically the marginal pdf of 
% - the first entry of the random vector Mu 
% - the (1,1)-entry of the random matrix inv(Sigma)
% when Mu and Sigma are jointly normal-inverse-Wishart: Mu ~ St(Mu_0,Sigma/T_0)
%                                                       inv(Sigma) ~ W(Nu_0,inv(Sigma_0)/Nu_0)
% See Ch.7 in "Risk and Asset Allocation" - Springer (2005), by A. Meucci

NumSimulations=size(Mu_Simul,1);
NumBins=round(10*log(NumSimulations));

figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mu

% plot empirical pdf (histogram)
subplot(2,1,1)
[n,x]=hist(Mu_Simul(:,1),NumBins);
D=x(2)-x(1);
n=n/(D*NumSimulations);
bar(x,n,1);

% superimpose analytical expectation
hold on 
h=plot(Mu_0(1),0,'.');
set(h,'color','r','markersize',15)

% superimpose analytical pdf
hold on 
x_lo=min(Mu_Simul(:,1));
x_hi=max(Mu_Simul(:,1));
x_grid=[x_lo: (x_hi-x_lo)/100 :x_hi];
m=Mu_0(1);
s=sqrt(Sigma_0(1,1)/T_0);
f=1/s*tpdf((x_grid-m)/s,Nu_0);
h=plot(x_grid,f,'r');

title([Legend '  Mu'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sigma

% plot empirical pdf (histogram)
subplot(2,1,2)
[n,x]=hist(InvSigma_Simul(:,1),NumBins);
D=x(2)-x(1);
n=n/(D*NumSimulations);
bar(x,n,1);

% superimpose analytical expectation
InvSigma_0=inv(Sigma_0);
hold on 
h=plot(InvSigma_0(1,1),0,'.');
set(h,'color','r','markersize',15)

% superimpose analytical pdf
hold on
x_lo=min(InvSigma_Simul(:,1));
x_hi=max(InvSigma_Simul(:,1));
x_grid=[x_lo: (x_hi-x_lo)/100 :x_hi];
sigma_square=InvSigma_0(1,1)/Nu_0;
A=Nu_0/2;
B=2*sigma_square;
f = gampdf(x_grid,A,B);
h=plot(x_grid,f,'r');

title([Legend '  inv(Sigma)'])