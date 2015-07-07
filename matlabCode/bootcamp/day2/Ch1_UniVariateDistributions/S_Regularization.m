% this script displyays the regularized empirical distribution 
% obtained by convolution with the smoothed Dirac delta, 
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci
% formulas (1.121) and (2.241)


clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
T=5;
Smoothness=.04;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Observations=normrnd(0,1,T,1); % generate observations (e.g. normal)

m=mean(Observations);
s=std(Observations);

% set range 
Step=(max(Observations)-min(Observations))/10000;
x=[min(Observations)-1000*Step : Step : max(Observations)+1000*Step];

% compute regularized pdf and cdf 
f_reg=0; 
F_reg=0;
for t=1:T
  f_reg = f_reg + 1/T*normpdf(x,Observations(t),Smoothness);
  F_reg = F_reg + 1/T*normcdf(x,Observations(t),Smoothness);
end
% compute regularized quantile (=inverse of regularized cdf, it is enough to transpose the plot)
u_reg = F_reg;
Q_reg = x; 


% compute true empirical cdf 
F_true=0*x;
for j=1:length(x)
   AreLess=0+(Observations<=x(j)); % vector of logical values: =1 when condition in () is satisfied, else =0.
                        % Adding 0+ converts logical variables to numbers
   Count=sum(AreLess);
   F_true(j)=Count/T;
end
% compute true empirical quantile
u_true = [.001 : .001 : .999];
OrderStats=sort(Observations);
Index=ceil(u_true*T);
Q_true = OrderStats(Index);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
figure
h=plot(x,f_reg,'b');
hold on
h=plot(Observations,0*Observations,'o');
grid on
title('pdf')

figure
h=plot(x,F_reg,'b');
hold on 
h=plot(x,F_true,'r');
grid on
legend('regularized', 'true','location','northwest')
title('cdf')

figure
h=plot(u_reg,Q_reg,'b');
hold on 
h=plot(u_true,Q_true,'r');
grid on
legend('regularized', 'true','location','northwest')
title('quantile')