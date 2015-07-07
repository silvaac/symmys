function PlotFrontier(Portfolios, vol)
% This function plots the efficient frontier in the plane of portfolio
% weights versus standard deviation. 

colors = {'c','w','m','g','b','r'};
numcolors = length(colors);

figure
[xx,N]=size(Portfolios);
Data=cumsum(Portfolios,2);
for n=1:N
    x=[1; [1 : xx]'; xx];
    x = [min(vol); vol; max(vol)];
    v = [min(vol); vol; max(vol)]*100;
    y=[0; Data(:,N-n+1); 0];
    hold on
    h=fill(v,y,colors{mod(n,numcolors)+1});
end
set(gca,'xlim',[v(1) v(end)],'ylim',[0 max(max(Data))])
xlabel('Risk %')
ylabel('Portfolio weights')
title('Efficient Frontier')
