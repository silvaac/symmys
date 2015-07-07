function [x_Hor,f_Hor,F_Hor]=ProjectionT(nu,m,s,T)
% this function performs the horizon projection of a t-distributed invariant
% inputs: 
% nu = t-distribution's degree of freedom
% s = t-distribution's scatter parameter
% m = t-distribution's location parameter
% T = multiple of the estimation period to the invesment horizon 
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci
% a detailed proof of the steps below can be found 
% at www.symmys.com > Book > Dowloads > Technical Appendices > Ch 3.7 "Numerical Market Projection"


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up grid
N=2^14; % coarseness level
a=-norminv(10^(-15),0,sqrt(T));
h=2*a/N;
Xi=[-a+h : h : a]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% discretized initial pdf (standardized)
f = 1/h*(tcdf(Xi+h/2,nu)-tcdf(Xi-h/2,nu));
f(N) = 1/h*(tcdf(-a+h/2,nu)-tcdf(-a,nu) + tcdf(a,nu)-tcdf(a-h/2,nu));

% discretized initial pdf (non-standardized)
X_Start=m+s*Xi;
f_Start=f/s;

% discretized characteristic function
Phi=fft(f);                     

% projection of discretized characteristic function
Signs=(-1).^([0:N-1]'*(T-1));   
Phi_T=h^(T-1)*Signs.*(Phi.^T);

% horizon discretized pdf (standardized)
f_T=ifft(Phi_T);

% horizon discretized pdf and cdf (non-standardized)
x_Hor=m*T+s*Xi;
f_Hor=f_T/s;
F_Hor=h*cumsum(f_Hor*s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots

% figure
% 
% subplot(3,1,1)
% u=plot(X_Start,f_Start);
% grid on
% title('estimation pdf')
% 
% subplot(3,1,2)
% u=plot(x_Hor,f_Hor);
% grid on
% title('horizon pdf')
% 
% subplot(3,1,3)
% u=plot(x_Hor,F_Hor);
% grid on
% title('horizon cdf')
% 
