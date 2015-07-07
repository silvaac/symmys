% this script generates random correlation matrices by randomizing the
% parameters of the general parametric specifications for a correlation matrix
% see "Risk and Asset Allocation"- Springer (2005), by A. Meucci

clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=3; % dimensionality of the problem
K=N*(N-1)/2;

J=10000; % number of simulations

% parametrization 1: randomize angles
Thetas=2*pi*rand(J,K);
Cosines=cos(Thetas);
Sines=sin(Thetas);

% parametrization 2: randomize sines/cosines
%Cosines=2*rand(J,K)-1;
%Sines=sqrt(1-Cosines.*Cosines);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute correlations in all scenarios
x=[0];
for n=2:N
    x=[x
        (n-1)*(n-2)/2 + 1];
end

CorrsAsTensor=zeros(J,N,N);
Eigs=[];
Z=zeros(N,1);
for j=1:J
    c=Cosines(j,:)';
    s=Sines(j,:)';

    V=Z;
    V(1)=1;
    for n=2:N
        v=Z;
        v(1)=c(x(n));
        S=[1; cumprod(s([x(n):x(n)+n-2]))];
        for k=2:n-1
            v(k)=S(k)*c(x(n)+k-1);
        end
        v(n)=S(end);

        V=[V v];
    end
    Corr=V'*V;
    CorrsAsTensor(j,:,:)=Corr;
    %Eigs=[Eigs
    %     eig(Corr)'];

end

% reassemble results in an entry-wise structure that runs on the upper
% triangular portion of the correlation

CorrsAsEntries=[];
k=0;
for n=1:N
    for m=n+1:N
        k=k+1;
        CorrsAsEntries(k).Values = CorrsAsTensor(:,n,m);
        CorrsAsEntries(k).Names = [num2str(n) '/' num2str(m)];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
Plots=1;  % choose to plot or not
if Plots

    % univariate marginals
    K=length(CorrsAsEntries);
    Nbins=round(20*log(J));
    for k=1:K
        figure
        hist(CorrsAsEntries(k).Values,Nbins);
        title(CorrsAsEntries(k).Names);
    end

    % bivariate marginals
    for k=1:K
        for j=k+1:K
            figure
            plot(CorrsAsEntries(k).Values,CorrsAsEntries(j).Values,'.');
            title([CorrsAsEntries(k).Names ' - ' CorrsAsEntries(j).Names]);
        end
    end
end