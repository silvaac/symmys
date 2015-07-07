clc; clear; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J=10000;
N=10;
sigmasq=1;

U=rand(J,1);
X=[];
for n=1:N
    a=n/2;
    b=2*sigmasq;
    X = [X gaminv(U,a,b)];
end