clear; clc; close all
% code to project a distribution in the future according to the i.i.d.-implied square-root rule
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci


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

Col=[.8 .8 .8];
subplot('Position',[.05 .55 .9 .4])
plot([0 ts],[0*ps' C.q],'r')
for i=1:length(ts)
    hold on
    xx=[ts(i) ; ts(i)+C.pdf(:,i) ; ts(i)];
    yy=[min(C.x(:,i)) ; C.x(:,i) ; max(C.x(:,i))];
    fill(xx,yy,Col)
end
xlim([0 max(xx)*1.01])
Ylim=[min(yy)*1.5 max(yy)*1.5];
ylim(Ylim)
grid on

subplot('Position',[.05 .05 .9 .40])
plot([0 ts],[0*ps' R.q],'r')
for i=1:length(ts)
    hold on
    xx=[ts(i) ; ts(i)+R.pdf(:,i) ; ts(i)];
    yy=[min(R.x(:,i)) ; R.x(:,i) ; max(R.x(:,i))];
    fill(xx,yy,Col)
end
xlim([0 max(xx)*1.01])
ylim(Ylim)
grid on
