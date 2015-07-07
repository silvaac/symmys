% this script generates a uniform sample on the interior of an ellipsoid in N dimensions
% See "Risk and Asset Allocation" - Springer (2005), by A. Meucci
% this script is intuitive but inefficient, see S_UniformEllipsoid2 for an efficient approach

clc; clear; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=3;
Mu=rand(N,1);
A=rand(N,N)-.5;
Sigma=A*A';

NumSimul=10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate uniform sample on the cube
VolSphere=pi^(N/2)/gamma(N/2+1);
VolCube=2^N;
NumSimulCube=ceil(NumSimul*VolCube/VolSphere * 1.1); % set # of simulations according to (A.78)
U=[];
for n=1:N
    u=2*rand(NumSimulCube,1)-1;
    U=[U u];
end

% keep ony observations within the sphere
Zsquare = sum(U.*U,2);
Keep=find(Zsquare<=1);
Y=U(Keep,:);
Y=Y(1:NumSimul,:);

% affine transformation from sphere to ellipsoid
X = ones(NumSimul,1)*Mu'+ Y*A';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scatter-plot the first two components
figure
plot(X(:,1),X(:,2),'.')
grid on
axis equal
