% this script shows the diversification effect of correlation in a simplified market of two variables
% see "Risk and Asset Allocation" - Springer (2005), by A. Meucci

clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_Low=0.05;
E_High=0.10;
S_Low=0.10;
S_High=0.25;
r=[-1  0  1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

% portfolio weights
w=[-0.5: 1/200 : 1.7];
for i=1:length(r)
    S=diag([S_Low S_High])*[1 r(i);r(i) 1]*diag([S_Low S_High]);
    m=[E_Low E_High]';
    E_r=[];
    Vol_r=[];
    for j= 1:length(w)
        u=[w(j) 1-w(j)]';
        E_r=[E_r
            u'*m];
        Vol_r=[Vol_r
            sqrt(u'*S*u)];
    end
    plot(Vol_r, E_r);
    hold on
end
grid on
