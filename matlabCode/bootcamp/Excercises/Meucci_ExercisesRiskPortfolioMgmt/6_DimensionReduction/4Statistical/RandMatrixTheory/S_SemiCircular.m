clc; close all; clear;

N=1000;  % matrix size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% empirical

%X=randn(N);                % normal
%X=(rand(N)-.5)*sqrt(12);   % uniform
X=log(rand(N))+1;           % exponential

Y=(X+X')/(2*sqrt(2*N));    % symmetrize and rescale
E=eig(Y)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% theoretical
t=[-1:.01:1]; 
g=2/pi*sqrt(1-t.^2);

NumBins=ceil(10*log(length(E)));
[b t_]=hist(E,NumBins);
D=t_(2)-t_(1);
h=b/(D*N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
bar(t_,h);
hold on;
h=plot(t,g,'r');
set(h,'linewidth',3)