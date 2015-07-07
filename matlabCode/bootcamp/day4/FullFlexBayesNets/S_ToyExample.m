% This toy example illustrates the use of Entropy Pooling to compute Fully Flexible Bayesian networks, see 
% A. Meucci (2010) "Fully Flexible Bayesian Networks", working paper
% 
%  Most recent version of article and code available at
%  http://www.symmys.com/node/152

clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up scenarios table and prior distribution

x_1=[1 2 3]; % {L,M,H}
x_2=[1 2];   % {D,S}
x_3=[1 2];   % {C,R}

X=[];
for i=1:length(x_1)
    for k=1:length(x_2)
        for l=1:length(x_3)
            X=[X
            x_1(i) x_2(k) x_3(l)];
        end
    end
end

[J,N]=size(X);
p=ones(J,1)/J;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input views
% statement: View(k).Who (e.g. [1 3])= View(k).Equal (e.g. {[2 3] [1 3 5]})
% optional conditional statement: View(k).Cond_Who (e.g. [2])= View(k).Cond_Equal (e.g. {[1]})
% amount of stress is quantified as Prob(statement) <= View(k).v if View(k).sgn = 1;
%                                   Prob(statement) >= View(k).v if View(k).sgn = -1;
% confidence in stress is quantified in View(k).c in (0,1)

View(1).Who=[1];  
View(1).Equal={[2 3]};
View(1).Cond_Who=[2]; 
View(1).Cond_Equal={[1]};
View(1).v = .7;
View(1).sgn =-1;
View(1).c=.5;

View(2).Who=[2];  
View(2).Equal={[1]};
View(2).Cond_Who=[]; 
View(2).Cond_Equal={[]};
View(2).v =.3;
View(2).sgn =-1;
View(2).c=.5;

% create linear constraint representation of views
[A,b,g] = CondProbViews(View,X);

% enforce consistence
db = Tweak(A,b,g);
b = b + db;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute posterior
Aeq = ones(1,J);  % constrain probabilities to sum to one
beq = 1;

% ...compute posterior probabilities
p_ = EntropyProg(p,A,b,Aeq ,beq);