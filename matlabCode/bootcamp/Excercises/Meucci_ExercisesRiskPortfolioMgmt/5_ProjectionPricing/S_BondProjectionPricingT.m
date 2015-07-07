% this script projects the distribution of the market invariants for the bond market 
% (i.e. the changes in yield to maturity) from the estimation interval (Student t assumption) 
% to the investment horizon using the FFT
% Then it computes the distribution of prices at the investment horizon by full Monte Carlo
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
tau=4/52;        % time to horizon expressed in years
tau_tilde=1/52;  % estimation period expressed in years

FlatCurve=.04;   
TimesToMat=[4 5 10 52 520]'/52; % time to maturity of selected bonds expressed in years

% determine the parameters of the distribution of the invariants (changes in yield to maturity)
Periods=tau/tau_tilde; % number of estimation periods until the investment horizon
u_minus_tau=TimesToMat-tau;

nu=8;
mus=0*u_minus_tau;
sigmas=(20+5/4*u_minus_tau)/10000;

Num_Scenarios=100000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projection and pricing 
BondCurrent_Prices_Shifted=exp(-FlatCurve*u_minus_tau);
BondCurrent_Prices=exp(-FlatCurve*TimesToMat);

% generate common source of randomness
U=rand(Num_Scenarios,1);  

N=length(TimesToMat); % number of bonds
for n=1:N
    % project bond market to horizon
    [x,f,F]=ProjectionT(nu,mus(n),sigmas(n),Periods);
    
    % generate co-dependent changes in yield-to-maturity
    DY_Scenarios = interp1(F,x,U,'linear','extrap'); 

    % compute the horizon prices, (3.81) in "Risk and Asset Allocation" - Springer
    X=-u_minus_tau(n)*DY_Scenarios;
    Z=BondCurrent_Prices_Shifted(n)*exp(X); 
    
    % compute and plot linear returns
    L=Z/BondCurrent_Prices(n)-1;  
    subplot(N,1,n)
    hist(L,  round(10*log(Num_Scenarios))  );
    xlabel(['Linear returns for bond ' num2str(n)])
end