% this script shows how the goodness of the Taylor approximation 
% (theta/delta/gamma - carry/duration/convexity) 
% for a pricing function depends on the natural scale of the risk factor
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; close all; clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
NumSimul=100000;

% specify the distribution of the invariants (can be any distribution, normal here)
X=.5+.6*trnd(20,NumSimul,1);

% set the order in the Taylor expansion
Order=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% exact pricing
P=sin(X);

% approximate pricing
m=mean(X); % Taylor's pivot: can pick an arbitrary value, 
           % but choices other than the mean will require higher orders
P_Approx=sin(m);
for k=1:Order
    switch mod(k,4)
        case 0 
            d=sin(m);
        case 1
            d=cos(m);
        case 2
            d=-sin(m);
        case 3
            d=-cos(m);
    end
    P_Approx = P_Approx + d/prod(1:k)*((X-m).^k);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumBins=round(10*log(NumSimul));

figure

% approximate pricing
subplot('Position',[.25 .3 .15 .6]) 
[n,y]=hist(P_Approx,NumBins);
n=n/(NumSimul*(y(2)-y(1)));
h2=barh(y,n);
set(h2,'facecolor','g','edgecolor','g')
[y_lim]=get(gca,'ylim');
set(gca,'xtick',[])
grid on
title(['order ' num2str(Order) ' approx.'])

% true pricing
subplot('Position',[.05 .3 .15 .6]) 
[n,y]=hist(P,NumBins);
n=n/(NumSimul*(y(2)-y(1)));
h1=barh(y,n);
set(h1,'facecolor','r','edgecolor','r')
set(gca,'ylim',y_lim)
set(gca,'xtick',[])
grid on
title(['exact'])

% underlying market
subplot('Position',[.45 .05 .45 .2]) 
[n,x]=hist(X,NumBins);
n=n/(NumSimul*(x(2)-x(1)));
h3=bar(x,n);
set(h3,'facecolor','b','edgecolor','b')
set(gca,'ytick',[])
[x_lim]=get(gca,'xlim');
grid on

% approximate pricing function
subplot('Position',[.45 .3 .45 .6]) 
plot(x,sin(x),'r')
hold on
y=interp1(X,P_Approx,x,'linear','extrap');
plot(x,y,'g')
set(gca,'xlim',x_lim,'ylim',y_lim)
grid on
legend('exact',['order ' num2str(Order) ' approx.'],'location','southwest')