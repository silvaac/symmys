% this script projects the distribution of the market invariants for the stock market (i.e. the compounded returns)
% from the estimation interval (normal assumption) to the investment horizon
% Then it computes the distribution of prices at the investment horizon
% analytically, by full Monte Carlo, and by delta/duration approximation
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
tau_tilde=1/52;  % estimation period expressed in years
sig=.4;
P_T=1;

NumScenarios=1000000;

taus=[1/252 1/52 1/12 1 2];    % times to horizon expressed in years
tauName = {'1 day','1 week','1 month','1 year','2 years'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:length(taus)  % loop projection/pricing over different times to horizon
    tau=taus(k)

    % exact simulation of horizon prices
    C_Ttau=normrnd(0,sqrt(sig*sig*tau),NumScenarios,1);
    P_Ttau=P_T*exp(C_Ttau);

    % compute analytical pdf
    p_lo=min(P_Ttau);
    p_hi=max(P_Ttau);
    p=[p_lo : (p_hi-p_lo)/1000 : p_hi];
    m=log(P_T);
    s=sqrt(sig*sig*tau);
    f=lognpdf(p,m,s);

    % compute approximate pdf
    f_approx=normpdf(p,P_T,sqrt(P_T*P_T*sig*sig*tau));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plots
    figure
    NumBins=round(10*log(NumScenarios));
    [count_bins,p_bins]=hist(P_Ttau,NumBins);
    Scale=(p_bins(2)-p_bins(1))*NumScenarios;
    f_bins=count_bins/Scale;
    bar(p_bins,f_bins)

    hold on
    plot(p,f,'r')

    hold on
    plot(p,f_approx,'g')
    xlabel('Price at the horizon')
    ylabel('pdf')
    title(['Time to horizon \tau = ' tauName{k}])
    
    legend('full Monte Carlo','analytical','delta/duration')
end