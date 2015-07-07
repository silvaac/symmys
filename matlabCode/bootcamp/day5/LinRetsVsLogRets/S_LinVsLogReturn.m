% SYMMYS - Last version of code and article available at http://symmys.com/node/141

% Differences between compounded (or logarithmic) and linear return 
% see Meucci, A. (2010) "Linear versus Compounded Returns:Common Pitfalls in Risk and Portfolio Management"
% GARP Risk Professional, April, pp. 49-51 


clear; clc; close all

% in general R=exp(C)-1. Furtheremore, here we assume C~N(m*t,s^2*t) 
m=.05;
s=.25;

ts=[.1 : .3 : 3];
ps=[.01 .99];

D=.7*min(abs(diff(ts)));
C.q=[];
C.x=[];
C.pdf=[];
R.q=[];
R.x=[];
R.pdf=[];
Steps=100;
for i=1:length(ts)
    t=ts(i);
    
    q=norminv(ps,m*t,s*sqrt(t))';
    
    x=min(q) : (max(q)-min(q))/Steps : max(q);
    pdf=normpdf(x,m*t,s*sqrt(t));
    pdf=pdf/max(pdf)*D;
    
    C.q=[C.q q];
    C.pdf=[C.pdf pdf'];
    C.x=[C.x x'];

    q=exp(q)-1;
    
    x=min(q) : (max(q)-min(q))/Steps : max(q);
    pdf=lognpdf(x+1,m*t,s*sqrt(t));
    pdf=pdf/max(pdf)*D;
    
    R.pdf=[R.pdf pdf'];
    R.x=[R.x x'];

 end
R.q=exp(C.q)-1;