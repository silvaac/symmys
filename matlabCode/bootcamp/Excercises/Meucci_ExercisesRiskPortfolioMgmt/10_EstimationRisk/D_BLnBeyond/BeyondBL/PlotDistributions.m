function PlotDistributions(X,p,p_)

[J,N]=size(X);
NBins=round(10*log(J));
for n=1:N

    figure 
    
    % set ranges
    xl=min(X(:,n));
    xh=max(X(:,n));
    
    % prior numerical
    subplot(2,1,1)
    h2=pHist(X(:,n),p,NBins);
    xlim([xl xh])
    y1 = ylim;
    title('prior')
    
    % posterior numerical
    subplot(2,1,2)
    h3=pHist(X(:,n),p_,NBins);
    xlim([xl xh])
    y2 = ylim;
    ylim([min(y1(1),y2(1)),max(y1(2),y2(2))]);
    title('posterior')
    
    subplot(2,1,1)
    ylim([min(y1(1),y2(1)),max(y1(2),y2(2))]);

end