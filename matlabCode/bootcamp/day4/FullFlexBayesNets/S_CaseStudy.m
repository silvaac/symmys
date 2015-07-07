% This case study uses Entropy Pooling to compute Fully Flexible Bayesian networks for risk management, see 
% A. Meucci (2010) "Fully Flexible Bayesian Networks", working paper
% 
%  Most recent version of article and code available at
%  http://www.symmys.com/node/152

clear; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% upload scenarios table and prior distribution for changes in
% SWAP2YR SWAP10YR CDXIG S&P500 DollarIndex Crude Gold VIX 10YRInflationSwapRate

load FreaqEst 

[J,N]=size(X); 
e=.01;
p=(1-e)*p+e*ones(J,1)/J;
[m,s,C]=ComputeMoments(X,p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input views
% statement: View(k).Who (e.g. [1 3])= View(k).Equal (e.g. {[2 3] [1 3 5]})
% optional conditional statement: View(k).Cond_Who (e.g. [2])= View(k).Cond_Equal (e.g. {[1]})
% amount of stress is quantified as Prob(statement) <= View(k).v if View(k).sgn = 1;
%                                   Prob(statement) >= View(k).v if View(k).sgn = -1;
% confidence in stress is quantified in View(k).c in (0,1)

for n=1:N
    k=2*n-1;
    View(k).Who=[n];  
    View(k).Equal={[-1]};
    View(k).Cond_Who=[]; 
    View(k).Cond_Equal={[]};
    View(k).v = .4;
    View(k).sgn =-1;
    View(k).c=.5;

    k=2*n;
    View(k).Who=[n];  
    View(k).Equal={[1]};
    View(k).Cond_Who=[]; 
    View(k).Cond_Equal={[]};
    View(k).v = .4;
    View(k).sgn =-1;
    View(k).c=.5;

end

k=2*N+1;
View(k).Who=[1 2 8];
View(k).Equal={[-1 0] [-1 0] [1]};
View(k).Cond_Who=[4]; 
View(k).Cond_Equal={[-1]};
View(k).v = .9;
View(k).sgn =-1;
View(k).c=.1;

% create linear constraint representation of views on probabilities
[A,b,g] = CondProbViews(View,X);

% add constraint for view on correlation
C_12_=.6;
New_A=(X(:,1).*X(:,2))';
New_b=s(1)*s(2)*C_12_ + m(1)*m(2);
New_g=-log(1-.1);

A=[A
    New_A];
b=[b
    New_b];
g=[g
    New_g];
    
% enforce consistency
db = Tweak(A,b,g);
b = b+db;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute posterior
Aeq = ones(1,J);  % constrain probabilities to sum to one
beq = 1;
p_ = EntropyProg(p,A,b,Aeq ,beq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots

figure
subplot(2,1,1)
bar(p)
YLim=get(gca,'ylim');
subplot(2,1,2)
bar(p_)
set(gca,'ylim',YLim);