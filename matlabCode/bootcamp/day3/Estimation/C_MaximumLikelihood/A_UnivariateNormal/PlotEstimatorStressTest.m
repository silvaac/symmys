function a2=PlotEstimatorStressTest(Losses,Inef2s,Bias2s,Error2s,P,PName,EvalName);

NumSimulations=size(Losses,1);
Steps=size(Losses,2);

figure

% plot loss
subplot(2,1,1)

% set axes to define limits
a1=gca;
Pos=get(a1,'Position');
D=(P(end)-P(1))/Steps;
X_Lim = [P(1)-D  P(end)+D];
set(a1,'xlim',X_Lim,'ytick',[])
xlabel(PName)
title(['loss for estimator of ' EvalName],'fontweight','bold');

Num_Bins=round(10*log(NumSimulations));
Gray=.5;
for s=1:Steps
    % superimpose rotated histograms of loss
    a=axes;
    Z=Losses(:,s);
    [n,xout]=hist(Z,Num_Bins);
    h=barh(xout,n,1);
    set(h,'FaceColor',Gray*[1 1 1],'EdgeColor',Gray*[1 1 1])

    if s==1 % set vertical limit
        Top=prctile(Z,95);
        Rescale=Top/max(Z);
        YLim=get(a,'ylim');
        YLim=YLim*Rescale;
    end
    set(a,'ylim',YLim);
    Pos_x= Pos(1)+(P(s)-X_Lim(1))/(X_Lim(2)-X_Lim(1))*Pos(3);
    set(a,'Position',[Pos_x  Pos(2)  1/(2*Steps) Pos(4)])
    axis off
end

% plot summary evaluators
subplot(2,1,2);
h=bar(P,Inef2s+Bias2s,'r');
hold on
h=bar(P,Inef2s,'w');
hold on
h=plot(P,Error2s,'.');
set(h,'color','k','markersize',10)
a2=gca;
set(a2,'xlim',X_Lim)
legend('sq.bias','sq.ineff.','sq.err.','location','best')
xlabel(PName)
title(['summary evaluation for estimator of ' EvalName],'fontweight','bold');
