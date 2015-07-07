% this script shows how the empirical distribution converges to the tre
% distribution as the number of simulations/observations increases (Glivenko-Cantelli)
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci
% see also the technical appendices at symmys.com>book>downloads


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clean up workspace
clear; close all;  clc; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate i.i.d. data (e.g. normal)
T=1000; % length of sample
mu=0;
sigma=1;
X=normrnd(mu,sigma,T,1);
figure  
h=plot(X,0*X,'d'); % plot realized values
set(h,'MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute empirical cdf at realized values
xx=sort(X);     % (check that setting T in the denominator yields 
FF=[1:T]/(T);   % the built-in function "ecdf")
% plot
hold on 
h=plot(xx,FF,'.'); % plot empirical cdf

% compute empirical cdf at pre-specified values
Refine_Range=1000;
m=min(X);
M=max(X);
Delta=(M-m)/(Refine_Range-1);
x=[m-.1*Refine_Range*Delta : Delta : M+.1*Refine_Range*Delta]; % pre-specified values
F=0*x;
for j=1:length(x)
   AreLess=0+(X<=x(j)); % vector of logical values: =1 when condition in () is satisfied, else =0.
                        % Adding 0+ converts logical variables to numbers
   Count=sum(AreLess);
   F(j)=Count/T;
end
% plot
hold on   % plot empirical cdf at realized values
h=plot(x,F);

% compute true limit cdf at realized values
F_limit=normcdf(x,mu,sigma);
% plot
hold on 
h=plot(x,F_limit,'r');
legend(h,'limit cdf','location','northwest')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot histogram
figure
NumBins=round(5*log(T));
[n,Ds]=hist(X,NumBins); % compute bin width
hist(X,NumBins) % plot histogram

% plot (rescaled) limit pdf
f_limit=normpdf(x,mu,sigma);
D=Ds(2)-Ds(1);
hold on
h=plot(x,T*D*f_limit,'r');
grid on
legend(h,'limit pdf','location','northwest')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% quantile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute empirical quantile at pre-specified values 
u=[.01 : .01 : .99];  % range of quantiles (values between zero and one)
Q=prctile(X,u*100);

% plot
figure 
plot(u,Q);
grid on
% plot inverse (tranpose) of empirical cdf
hold on 
plot(FF,xx,'.');
grid on

% compute true limit quantile 
Q_limit=norminv(u,mu,sigma);
% plot
hold on 
h=plot(u,Q_limit,'r');
legend(h,'limit quantile','location','northwest')
