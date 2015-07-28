clear all
close all
load  'CaseStudy1.mat';

b_grid = [-1:0.5:1];
min_X= 100*min(X_path(1:tau_-1));
Delta_X = (100*mu_LT-min_X)/2;
X_grid = [min_X:Delta_X:min_X+Delta_X*4];
max_X= 100*max(X_path(1:tau_-1));
Delta_X = (-100*mu_view+max_X)/2;
X_grid_viewMean = sort([max_X:-Delta_X:max_X-Delta_X*4]);

f1 = figure(1);
set(f1,'color','w','units','normalized','position',[0.2 0.1 0.55 0.7]);
h1 = subplot(2,2,1,'Parent',f1);
ax1_h1 = gca;
set(h1,'units','normalized');
set(h1,'position',get(h1,'position').*[0.95 1 1 1]);
axis(h1,[0 t(tau_-1) min(b_grid) max(b_grid)]);
set(h1,'units','normalized','YTick',b_grid,'YTickLabel', num2str(b_grid','%1.2f'),'fontsize',8,'FontWeight', 'bold');
ylabel(h1,'Exposure (x10^4)','fontsize',10,'FontWeight', 'bold', 'Units','pixels');     
line(t(1:tau_-1),b_NoMI_prior/10^4,'Color','k','Parent',ax1_h1);
set(h1,'XColor','k','YColor','k')
ax2_h1 = axes('Position',get(ax1_h1,'Position'),...
            'XAxisLocation','top',...    
            'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','r');
line(t(1:tau_-1), repmat(100*mu_LT,1,tau_-1),'Color','r','Linestyle','--','Parent',ax2_h1);
line(t(1:tau_-1),100*X_path(1:tau_-1),'Color','r','Parent',ax2_h1);       
axis(ax2_h1,[0 t(tau_-1) min(X_grid) max(X_grid)]);
set(ax2_h1,'units','normalized','YTick',X_grid,'YTickLabel', num2str(X_grid','%1.2f'),'XTickLabel', [],'fontsize',8,'FontWeight', 'bold');
ylabel(ax2_h1,'10y rate (%)','fontsize',10,'FontWeight', 'bold', 'Units','pixels');     
title('$b_t^{Prior}$','interpreter','latex','fontsize',14,'FontWeight', 'bold')
legend('Long term mean','Location','NOrthWest')
%%%
h2 = subplot(2,2,2,'Parent',f1);
ax1_h2 = gca;
set(h2,'units','normalized');
set(h2,'position',get(h2,'position').*[1.05 1 1 1]);
line(t(1:tau_-1),b_NoMI_post/10^4,'Color','k','Parent',ax1_h2);
line(t(1:tau_-1),b_MI_Bellman_post(1:end-1,1)/10^4,'Color','b','Parent',ax1_h2);
axis(h2,[0 t(tau_-1) min(b_grid) max(b_grid)]);
set(h2,'units','normalized','YTick',b_grid,'YTickLabel', num2str(b_grid','%1.2f'),'fontsize',8,'FontWeight', 'bold');
ylabel(h2,'Exposure (x10^4)','fontsize',10,'FontWeight', 'bold', 'Units','pixels');     
set(h2,'XColor','k','YColor','k')
legend('No MI', 'With MI', 'location','SouthEast')
title('$b_t^{\ast}$','interpreter','latex','fontsize',14,'FontWeight', 'bold')
box on
%%%
h3 = subplot(2,2,3,'Parent',f1);
ax1_h3 = gca;
set(h3,'units','normalized');
set(h3,'position',get(h3,'position').*[0.95 1 1 1]);
axis(h3,[0 t(tau_-1) min(b_grid) max(b_grid)]);
set(h3,'units','normalized','YTick',b_grid,'YTickLabel', num2str(b_grid','%1.2f'),'fontsize',8,'FontWeight', 'bold');
ylabel(h3,'Exposure (x10^4)','fontsize',10,'FontWeight', 'bold', 'Units','pixels');     
xlabel('t(years)')
line(t(1:tau_-1),b_NoMI_LongTerm/10^4,'Color','k','Parent',ax1_h3);
set(h3,'XColor','k','YColor','k')
ax2_h3 = axes('Position',get(ax1_h3,'Position'),...
            'XAxisLocation','top',...    
            'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','r');
line(t(1:tau_-1), repmat(100*mu_LT,1,tau_-1),'Color','r','Linestyle','--','Parent',ax2_h3);
line(t(1:tau_-1),100*X_path(1:tau_-1),'Color','r','Parent',ax2_h3);       
axis(ax2_h3,[0 t(tau_-1) min(X_grid) max(X_grid)]);
set(ax2_h3,'units','normalized','YTick',X_grid,'YTickLabel', num2str(X_grid','%1.2f'),'XTickLabel', [],'fontsize',8,'FontWeight', 'bold');
ylabel(ax2_h3,'10y rate (%)','fontsize',10,'FontWeight', 'bold', 'Units','pixels');     
title('$b_t^{LongTerm}$','interpreter','latex','fontsize',14,'FontWeight', 'bold')
legend('Long term mean','location','NOrthWest')
%%%
h4 = subplot(2,2,4,'Parent',f1);
ax1_h4 = gca;
set(h4,'units','normalized');
set(h4,'position',get(h4,'position').*[1.05 1 1 1]);
axis(h4,[0 t(tau_-1) min(b_grid) max(b_grid)]);
set(h4,'units','normalized','YTick',b_grid,'YTickLabel', num2str(b_grid','%1.2f'),'fontsize',8,'FontWeight', 'bold');
ylabel(h4,'Exposure (x10^4)','fontsize',10,'FontWeight', 'bold', 'Units','pixels');     
line(t(1:tau_-1),b_NoMI_viewMean/10^4,'Color','k','Parent',ax1_h4);
set(h4,'XColor','k','YColor','k')
xlabel('t(years)')
ax2_h4 = axes('Position',get(ax1_h4,'Position'),...
            'XAxisLocation','top',...    
            'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','r');
line(t(1:tau_-1), repmat(100*mu_view,1,tau_-1),'Color','r','Linestyle','--','Parent',ax2_h4);
line(t(1:tau_-1),100*X_path(1:tau_-1),'Color','r','Parent',ax2_h4);       
axis(ax2_h4,[0 t(tau_-1) min(X_grid_viewMean) max(X_grid_viewMean)]);
set(ax2_h4,'units','normalized','YTick',X_grid_viewMean,'YTickLabel', num2str(X_grid_viewMean','%1.2f'),'XTickLabel', [],'fontsize',8,'FontWeight', 'bold');
ylabel(ax2_h4,'10y rate (%)','fontsize',10,'FontWeight', 'bold', 'Units','pixels');     
title('$b_t^{viewMean}$','interpreter','latex','fontsize',14,'FontWeight', 'bold')
legend('view mean level','location','SouthEast')

