close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% absolute coordinates: variance - expected value
figure
h=plot(TotRetCurve_Vs_1,TotRetCurve_Es_1,'r');
set(h,'linewidth',2);
hold on
h=plot(TotRetCurve_Vs_2,TotRetCurve_Es_2,'g');
set(h,'linewidth',2);
hold on
h=plot(TotRetCurve_Vs_3,TotRetCurve_Es_3,'k');
set(h,'linewidth',2);

hold on
h=plot(BenchRelCurve_Vs_1,BenchRelCurve_Es_1,'r');
set(h,'linewidth',2);
hold on
h=plot(BenchRelCurve_Vs_2,BenchRelCurve_Es_2,'g');
set(h,'linewidth',2);
hold on
h=plot(BenchRelCurve_Vs_3,BenchRelCurve_Es_3,'k');
set(h,'linewidth',2);

hold on
h=plot(MV_Variance,MV_ExpectedValue,'.','color','k');
set(h,'markersize',20);
hold on
h=plot(Sh_Variance,Sh_ExpectedValue,'.','color','k');
set(h,'markersize',20);
hold on
h=plot(Bench_Variance,Bench_ExpectedValue,'.','color','k');
set(h,'markersize',20);
xlim([0 max(BenchRelCurve_Vs_3)])
grid on

xlabel('variance')
ylabel('expected value')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% absolute coordinates: standard deviation - expected value
figure
h=plot(sqrt(TotRetCurve_Vs_1),TotRetCurve_Es_1,'r');
set(h,'linewidth',2);
hold on
h=plot(sqrt(TotRetCurve_Vs_2),TotRetCurve_Es_2,'g');
set(h,'linewidth',2);
hold on
h=plot(sqrt(TotRetCurve_Vs_3),TotRetCurve_Es_3,'k');
set(h,'linewidth',2);

hold on
h=plot(sqrt(BenchRelCurve_Vs_1),BenchRelCurve_Es_1,'r');
set(h,'linewidth',2);
hold on
h=plot(sqrt(BenchRelCurve_Vs_2),BenchRelCurve_Es_2,'g');
set(h,'linewidth',2);
hold on
h=plot(sqrt(BenchRelCurve_Vs_3),BenchRelCurve_Es_3,'k');
set(h,'linewidth',2);

hold on
h=plot(sqrt(MV_Variance),MV_ExpectedValue,'.','color','k');
set(h,'markersize',20);
hold on
h=plot(sqrt(Sh_Variance),Sh_ExpectedValue,'.','color','k');
set(h,'markersize',20);
hold on
h=plot(sqrt(Bench_Variance),Bench_ExpectedValue,'.','color','k');
set(h,'markersize',20);
xlim([0 sqrt(max(BenchRelCurve_Vs_3))])
grid on

xlabel('standard deviation')
ylabel('expected value')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% relative coordinates: squared tracking error - expected overperformance

figure % squared 
h=plot(TotRetCurve_TE2s_1,TotRetCurve_EOPs_1,'r');
set(h,'linewidth',2);
hold on
h=plot(TotRetCurve_TE2s_2,TotRetCurve_EOPs_2,'g');
set(h,'linewidth',2);
hold on
h=plot(TotRetCurve_TE2s_3,TotRetCurve_EOPs_3,'k');
set(h,'linewidth',2);

hold on
h=plot(BenchRelCurve_TE2s_1,BenchRelCurve_EOPs_1,'r');
set(h,'linewidth',2);
hold on
h=plot(BenchRelCurve_TE2s_2,BenchRelCurve_EOPs_2,'g');
set(h,'linewidth',2);
hold on

h=plot(BenchRelCurve_TE2s_3,BenchRelCurve_EOPs_3,'k');
set(h,'linewidth',2);

hold on
h=plot(MV_TE2,MV_EOP,'.','color','k');
set(h,'markersize',20);
hold on
h=plot(Sh_TE2,Sh_EOP,'.','color','k');
set(h,'markersize',20);
hold on
h=plot(Bench_TE2,Bench_EOP,'.','color','k');
set(h,'markersize',20);
grid on 

xlabel('squared tracking error')
ylabel('expected outperformance')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% relative coordinates: tracking error - expected overperformance
figure 
h=plot(sqrt(TotRetCurve_TE2s_1),TotRetCurve_EOPs_1,'r');
set(h,'linewidth',2);
hold on
h=plot(sqrt(TotRetCurve_TE2s_2),TotRetCurve_EOPs_2,'g');
set(h,'linewidth',2);
hold on
h=plot(sqrt(TotRetCurve_TE2s_3),TotRetCurve_EOPs_3,'k');
set(h,'linewidth',2);

hold on
h=plot(sqrt(BenchRelCurve_TE2s_1),BenchRelCurve_EOPs_1,'r');
set(h,'linewidth',2);
hold on
h=plot(sqrt(BenchRelCurve_TE2s_2),BenchRelCurve_EOPs_2,'g');
set(h,'linewidth',2);
hold on
h=plot(sqrt(BenchRelCurve_TE2s_3),BenchRelCurve_EOPs_3,'k');
set(h,'linewidth',2);

hold on
h=plot(sqrt(MV_TE2),MV_EOP,'.','color','k');
set(h,'markersize',20);
hold on
h=plot(sqrt(Sh_TE2),Sh_EOP,'.','color','k');
set(h,'markersize',20);
hold on
h=plot(sqrt(Bench_TE2),Bench_EOP,'.','color','k');
set(h,'markersize',20);
grid on 

xlabel('tracking error')
ylabel('expected outperformance')
