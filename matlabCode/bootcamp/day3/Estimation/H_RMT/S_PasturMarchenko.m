clc; close all; clear;

T=2000;
N=800;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% empirical

%X=(rand(T,N)-.5)*sqrt(12);
%X=randn(T,N);
%S=X'*X/T;
%kk=rand(N)-.5; kk=kk*kk'; [a,Sig]=cov2corr(kk);
Sig=eye(N);
S=wishrnd(Sig,T)/T;
E=eig(S);

NumBins=ceil(10*log(length(E)));
[B x]=hist(E,NumBins,1); 
D=x(2)-x(1);
p=B/(D*N);
bar(x,p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% theoretical
q=N/T;
l_min=(1-sqrt(q))^2;
l_max=(1+sqrt(q))^2;

l=[x(1) : (x(end)-x(1))/100 : x(end)];
a=max(l_max-l,0);
b=max(l-l_min,0);
y=1./(q.*2*pi*l).*sqrt(a.*b);
hold on
plot(l,y,'r');