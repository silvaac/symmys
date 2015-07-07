clc; close all; clear;

N=100;  % matrix size
J=1000; % samples

E=[]; % record eigenvalues
for j=1:J,
    A=randn(N);
    %A=(rand(N)-.5)*sqrt(12);
    A=(A+A')/(2*sqrt(2*N));  
  e=eig(A)'; 
  E=[E e];
end

NumBins=ceil(10*log(length(E)));
[B x]=hist(E,NumBins); 
D=x(2)-x(1);
p=B/(D*N*J);
bar(x,p);
hold on;
t=[-1:.01:1]; plot(t,2/pi*sqrt(1-t.^2),'r');
