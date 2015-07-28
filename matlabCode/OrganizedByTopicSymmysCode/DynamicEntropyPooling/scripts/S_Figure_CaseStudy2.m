clear all
close all
load  'CaseStudy1.mat';
clearvars -except b_MI_Bellman_post
b_MI_Bellman_post_1view = b_MI_Bellman_post;
clearvars b_MI_Bellman_post
load  'CaseStudy2.mat';

b_grid = [-1.2:0.6:1.2];

min_XLT= 100*min(X_path(1:tau_-1,1));
Delta_XLT = (100*mu_LT(1)-min_XLT)/2;
X_gridLongTerm = [min_XLT:Delta_XLT:min_XLT+Delta_XLT*4] ;

max_TIPLT= 100*max(TIP_path(1:tau_-1,1));
Delta_TIPLT = (-100*mu_LT(2) + max_TIPLT)/2;
TIP_gridLongTerm = sort([max_TIPLT:-Delta_TIPLT:max_TIPLT-Delta_TIPLT*4]);

min_TIPview= 100*min(TIP_path(1:tau_-1,1));
Delta_TIPview = (+100*mu_view(2) - min_TIPview)/2;
TIP_gridView = sort([min_TIPview:Delta_TIPview:min_TIPview+Delta_TIPview*4]);

max_Xview= 100*max(X_path(1:tau_-1,1));
Delta_Xview = (-100*mu_view(1) + max_Xview)/2;
X_gridView = sort([max_Xview:-Delta_Xview:max_Xview-Delta_Xview*4]);

close all
f1 = figure(1);
set(f1,'color','w','units','normalized','position',[0.2 0.1 0.55 0.8]);
h1 = subplot(3,2,1,'Parent',f1);
ax1_h1 = gca;
set(h1,'units','normalized');
set(h1,'position',get(h1,'position').*[0.95 1 1 1]);
axis(h1,[0 t(tau_-1) min(b_grid) max(b_grid)]);
set(h1,'units','normalized','YTick',b_grid,'YTickLabel', num2str(b_grid','%1.1f'),'fontsize',8,'FontWeight', 'bold');
ylabel(h1,'Exposure (x10^4)','fontsize',10,'FontWeight', 'bold', 'Units','pixels');     
line(t(1:tau_-1),b_NoMI_LongTermX/10^4,'Color','k','Parent',ax1_h1);
set(h1,'XColor','k','YColor','k')
ax2_h1 = axes('Position',get(ax1_h1,'Position'),...
            'XAxisLocation','top',...    
            'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','r');
axis(ax2_h1,[0 t(tau_-1) min(X_gridLongTerm) max(X_gridLongTerm)]);
set(ax2_h1,'units','normalized','YTick',X_gridLongTerm,'YTickLabel', num2str(X_gridLongTerm','%1.1f'),'XTickLabel', [],'fontsize',8,'FontWeight', 'bold');
line(t(1:tau_-1), repmat(100*mu_LT(1),1,tau_-1),'Color','r','Linestyle','--','Parent',ax2_h1);
line(t(1:tau_-1),100*X_path(1:tau_-1,1),'Color','r','Parent',ax2_h1);       
ylabel(ax2_h1,'10y rate (%)','fontsize',10,'FontWeight', 'bold', 'Units','pixels');     
title('$b_t^{LongTerm,X_1}$','interpreter','latex','fontsize',12,'FontWeight', 'bold')
legend('Long term level','location','NorthEast')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h2 = subplot(3,2,2,'Parent',f1);
ax1_h2 = gca;
set(h2,'units','normalized');
set(h2,'position',get(h2,'position').*[1.05 1 1 1]);
axis(h2,[0 t(tau_-1) min(b_grid) max(b_grid)]);
set(h2,'units','normalized','YTick',b_grid,'YTickLabel', num2str(b_grid','%1.1f'),'fontsize',8,'FontWeight', 'bold');
ylabel(h2,'Exposure (x10^4)','fontsize',10,'FontWeight', 'bold', 'Units','pixels');     
line(t(1:tau_-1),b_NoMI_LongTermTIP/10^4,'Color','k','Parent',ax1_h2);
set(h2,'XColor','k','YColor','k')
ax2_h2 = axes('Position',get(ax1_h2,'Position'),...
            'XAxisLocation','top',...    
            'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','g');
line(t(1:tau_-1), repmat(100*mu_LT(2),1,tau_-1),'Color','g','Linestyle','--','Parent',ax2_h2);
line(t(1:tau_-1),100*X_path(1:tau_-1,2),'Color','g','Parent',ax2_h2);       
axis(ax2_h2,[0 t(tau_-1) min(TIP_gridLongTerm) max(TIP_gridLongTerm)]);
set(ax2_h2,'units','normalized','YTick',TIP_gridLongTerm,'YTickLabel', num2str(TIP_gridLongTerm','%1.1f'),'XTickLabel', [],'fontsize',8,'FontWeight', 'bold');
ylabel(ax2_h2,'TIP spread (%)','fontsize',10,'FontWeight', 'bold', 'Units','pixels');     
title('$b_t^{LongTerm,X_2}$','interpreter','latex','fontsize',12,'FontWeight', 'bold')
legend('Long term level','location','SouthEast')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h3 = subplot(3,2,3,'Parent',f1);
ax1_h3 = gca;
set(h3,'units','normalized');
set(h3,'position',get(h3,'position').*[0.95 1 1 1]);
axis(h3,[0 t(tau_-1) min(b_grid) max(b_grid)]);
set(h3,'units','normalized','YTick',b_grid,'YTickLabel', num2str(b_grid','%1.1f'),'fontsize',8,'FontWeight', 'bold');
ylabel(h3,'Exposure (x10^4)','fontsize',10,'FontWeight', 'bold', 'Units','pixels');     
line(t(1:tau_-1),b_NoMI_viewMeanX/10^4,'Color','k','Parent',ax1_h3);
set(h3,'XColor','k','YColor','k')
%xlabel('t(years)')
ax2_h3 = axes('Position',get(ax1_h3,'Position'),...
            'XAxisLocation','top',...    
            'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','r');
line(t(1:tau_-1), repmat(100*mu_view(1),1,tau_-1),'Color','r','Linestyle','--','Parent',ax2_h3);
line(t(1:tau_-1),100*X_path(1:tau_-1,1),'Color','r','Parent',ax2_h3);       
axis(ax2_h3,[0 t(tau_-1) min(X_gridView) max(X_gridView)]);
ylabel(h3,'Exposure (x10^4)','fontsize',10,'FontWeight', 'bold', 'Units','pixels');     
set(ax2_h3,'units','normalized','YTick',X_gridView,'YTickLabel', num2str(X_gridView','%1.1f'),'XTickLabel', [],'fontsize',8,'FontWeight', 'bold');
ylabel(ax2_h3,'10y rate (%)','fontsize',10,'FontWeight', 'bold', 'Units','pixels');     
title('$b_t^{ViewMean,X_1}$','interpreter','latex','fontsize',12,'FontWeight', 'bold')
legend('View level','location','SouthEast')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h4 = subplot(3,2,4,'Parent',f1);
ax1_h4 = gca;
set(h4,'units','normalized');
set(h4,'position',get(h4,'position').*[1.05 1 1 1]);
axis(h4,[0 t(tau_-1) min(b_grid) max(b_grid)]);
set(h4,'units','normalized','YTick',b_grid,'YTickLabel', num2str(b_grid','%1.1f'),'fontsize',8,'FontWeight', 'bold');
ylabel(h4,'Exposure (x10^4)','fontsize',10,'FontWeight', 'bold', 'Units','pixels');     
line(t(1:tau_-1),b_NoMI_viewMeanTIP/10^4,'Color','k','Parent',ax1_h4);
set(h4,'XColor','k','YColor','k')
%xlabel('t(years)')
ax2_h4 = axes('Position',get(ax1_h4,'Position'),...
            'XAxisLocation','top',...    
            'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','g');
line(t(1:tau_-1), repmat(100*mu_view(2),1,tau_-1),'Color','g','Linestyle','--','Parent',ax2_h4);
line(t(1:tau_-1),100*X_path(1:tau_-1,2),'Color','g','Parent',ax2_h4);       
axis(ax2_h4,[0 t(tau_-1) min(TIP_gridView) max(TIP_gridView)]);
ylabel(h4,'Exposure (x10^4)','fontsize',10,'FontWeight', 'bold', 'Units','pixels');     
set(ax2_h4,'units','normalized','YTick',TIP_gridView,'YTickLabel', num2str(TIP_gridView','%1.1f'),'XTickLabel', [],'fontsize',8,'FontWeight', 'bold');
ylabel(ax2_h4,'TIP spread (%)','fontsize',10,'FontWeight', 'bold', 'Units','pixels');     
title('$b_t^{ViewMean,X_2}$','interpreter','latex','fontsize',12,'FontWeight', 'bold')
legend('view level','Location','SouthEast')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h5 = subplot(3,2,5,'Parent',f1);
ax1_h5 = gca;
set(h5,'units','normalized');
set(h5,'position',get(h5,'position').*[0.95 1 1 1]);
bpost_grid = b_grid;
axis(h5,[0 t(tau_-1) min(b_grid) max(b_grid)]);
set(h5,'units','normalized','YTick',bpost_grid,'YTickLabel', num2str(bpost_grid','%1.1f'),'fontsize',8,'FontWeight', 'bold');
line(t(1:tau_-1),b_NoMI_post(:,1)/10^4,'Color','k','Parent',ax1_h5);
line(t(1:tau_-1),b_MI_Bellman_post(1:end-1,1)/10^4,'Color','b','Parent',ax1_h5);
ylabel(h5,'Exposure (x10^4)','fontsize',10,'FontWeight', 'bold', 'Units','pixels');     
set(h5,'XColor','k','YColor','k')
xlabel('t(years)')
legend('No MI', 'With MI','location','SouthEast')
title('$b_t^{\ast}$','interpreter','latex','fontsize',12,'FontWeight', 'bold')
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h6 = subplot(3,2,6,'Parent',f1);
ax1_h6 = gca;
set(h6,'units','normalized');
set(h6,'position',get(h6,'position').*[1.05 1 1 1]);
bpost_grid = b_grid;
axis(h6,[0 t(tau_-1) min(b_grid) max(b_grid)]);
set(h6,'units','normalized','YTick',bpost_grid,'YTickLabel', num2str(bpost_grid','%1.1f'),'fontsize',8,'FontWeight', 'bold');
line(t(1:tau_-1),b_NoMI_prior(:,1)/10^4,'Color','k','LineStyle',':','Parent',ax1_h6);
line(t(1:tau_-1),b_MI_Bellman_post_1view(1:end-1,1)/10^4,'Color','k','Parent',ax1_h6);
line(t(1:tau_-1),b_MI_Bellman_post(1:end-1,1)/10^4,'Color','b','Parent',ax1_h6);
ylabel(h6,'Exposure (x10^4)','fontsize',10,'FontWeight', 'bold', 'Units','pixels');     
set(h6,'XColor','k','YColor','k')
xlabel('t(years)')
legend('prior', '1 view','2 views','location','SouthEast')
title('Comparison of $b_t^{\ast}$','interpreter','latex','fontsize',12,'FontWeight', 'bold')
box on

