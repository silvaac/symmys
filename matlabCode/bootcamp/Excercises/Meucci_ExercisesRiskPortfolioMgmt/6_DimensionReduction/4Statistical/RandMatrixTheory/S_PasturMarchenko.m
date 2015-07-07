clc; close all; clear;

T=1800;
N=700;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% empirical

%X=randn(T,N);                % normal
%X=(rand(T,N)-.5)*sqrt(12);   % uniform
X=log(rand(T,N))+1;           % exponential

Y=(X'*X)/T;    % symmetrize and rescale
E=eig(Y)';

NumBins=ceil(10*log(length(E)));
[b t_]=hist(E,NumBins,1); 
D=t_(2)-t_(1);
h=b/(D*N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% theoretical
q=N/T;
t_min=(1-sqrt(q))^2;
t_max=(1+sqrt(q))^2;
t=[t_(1) : (t_(end)-t_(1))/100 : t_(end)];
a=max(t_max-t,0);
b=max(t-t_min,0);
y=1./(q*2*pi*t).*sqrt(a.*b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
bar(t_,h);
hold on
h=plot(t,y,'r');
set(h,'linewidth',3)