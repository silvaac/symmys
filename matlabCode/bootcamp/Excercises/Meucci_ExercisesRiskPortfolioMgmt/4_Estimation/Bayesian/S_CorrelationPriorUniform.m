% this script shows how a jointly uniform prior on the correlations implies  
% that the marginal distribution of each correlation is peaked around zero.
% see "Risk and Asset Allocation"- Springer (2005), by A. Meucci

clear; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=3; % dimensionality of the problem
K=N*(N-1)/2;

J=10000; % number of simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute correlations in all scenarios

CorrsAsTensor=zeros(J,N,N);
Eigs=[];
j=1;
tic
while j<J
    C=2*rand(1,K)-1;
    Corr=eye(N);
    k=0;
    for n=1:N
        for m=n+1:N
            k=k+1;
            Corr(n,m)=C(k);
            Corr(m,n)=Corr(n,m);
        end
    end
    E=eig(Corr);
    if min(E)>0
        CorrsAsTensor(j,:,:)=Corr;
        j=j+1;
    end
    Eigs=[Eigs
        E'];

end
toc

% reassemble results in an entry-wise structure that runs on the upper
% triangular portion of the correlation
CorrsAsEntries=[];
k=0;
for n=1:N
    for m=n+1:N
        k=k+1;
        CorrsAsEntries(k).Values = CorrsAsTensor(:,n,m);
        CorrsAsEntries(k).Names = ['\theta_{' num2str(n) num2str(m) '}'];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots

% univariate marginals
K=length(CorrsAsEntries);
Nbins=round(5*log(J));
for k=1:K
    figure
    hist(CorrsAsEntries(k).Values,Nbins);
    title(sprintf('Histogram of %s',CorrsAsEntries(k).Names));
end

break
% bivariate marginals
for k=1:K
    for j=k+1:K
        figure
        plot(CorrsAsEntries(k).Values,CorrsAsEntries(j).Values,'.');
        title([CorrsAsEntries(k).Names ' - ' CorrsAsEntries(j).Names]);
    end
end
