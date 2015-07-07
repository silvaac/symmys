% this script decomposes the bivariate normal distribution into its radial and uniform components
% then it uses the uniform component to generate a spherically symmetric distribution
% see Sec. 2.7.1 in "Risk and Asset Allocation"-Springer (2005), by A. Meucci

clc; clear; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumSimul=10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = mvnrnd([0 0]',eye(2),NumSimul);
Norm_Y = sqrt(sum(Y.*Y,2));  

U = Y./(Norm_Y *[1 1]);       % uniform distribution on circle

% arbitrary radial distribution (choose one or set your own)
R_Arbitrary = chi2rnd(30,NumSimul,1); 
%R_Arbitrary = lognrnd(0,1,NumSimul,1); % arbitrary bimodal radial distribution
%R_Arbitrary = rand(NumSimul,1); % arbitrary bimodal radial distribution
X = (R_Arbitrary*[1 1]).*U;   % bi-variate spherically symmetric distribution: X = RU

% generate random rotations
for n=1:3
    t=2*pi*rand();
    E=[cos(t) -sin(t)
       sin(t)  cos(t)];
   Rotation(n).X=X*E';
   Rotation(n).t=t;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure % original
subplot(2,2,1)
plot(X(:,1),X(:,2),'.');
axis equal
grid on
xlabel('Original')
Xlim=get(gca,'xlim');
Ylim=get(gca,'ylim');

% rotations
for n=1:3
    subplot(2,2,n+1)
    h=plot(Rotation(n).X(:,1),Rotation(n).X(:,2),'.');
    set(gca,'xlim',Xlim','ylim',Ylim)
    grid on
    xlabel(['Rotation = ' num2str(Rotation(n).t/pi) ' pi'])
end
