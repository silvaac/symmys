% this script shows how the empirical distribution converges to the tre
% distribution as the number of simulations/observations increases (Glivenko-Cantelli)
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clean up workspace
clear; close all;  clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate i.i.d. data (e.g. normal)
T=100; % length of sample
Mu=[1 1];
Sigma=[1 .5; .5 1];
X=mvnrnd(Mu,Sigma,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



figure

% plot histogram
NBins =round(sqrt(T)/5);
[N,Ds]=hist3(X,[NBins NBins]); % compute bin width
hist3(X,[NBins NBins]) % plot histogram

% plot (rescaled) limit pdf
hold on
Ds_a=Ds{1}';  % x1-coordinates of bin centers
Ds_b=Ds{2}';  % x2-coordinates of bin centers
D_a=Ds_a(2)-Ds_a(1); % bin width on x1-coordinate
D_b=Ds_b(2)-Ds_b(1); % bin width on x2-coordinate
for i=1:length(Ds_a)
    for j=1:length(Ds_b)
        x_1=Ds_a(i);
        x_2=Ds_b(j);
        f_limit(i,j) = mvnpdf([x_1 x_2],Mu,Sigma);

    end
end
[X_1,X_2]=meshgrid(Ds_a,Ds_b);
mesh(X_1,X_2,T*D_a*D_b*f_limit);