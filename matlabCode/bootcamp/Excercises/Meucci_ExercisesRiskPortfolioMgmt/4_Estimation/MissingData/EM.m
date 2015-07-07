function [E_EM, S_EM, Y, CountLoop]=EM(X)
% this function implements the Expectation-Maximization (EM) algorithm to recover 
% missing observations in a time series 
% see "Risk and Asset Allocation"- Springer (2005), by A. Meucci

[T,N]=size(X);

% E-M initialization
X_Init=[];
for t=1:T
    if isempty(X(isnan(X(t,:))))
        X_Init=[X_Init
            X(t,:)];
    end
end
M=mean(X_Init)';
S=cov(X_Init);

Tolerance=10^(-6)*mean(  [M
                          sqrt(diag(S))]  );

% E-M loop
Convergence=0;
CountLoop=0;
Y=X;
while ~Convergence
    CountLoop=CountLoop+1;

    % Step 1: estimation
    C=zeros(T,N,N);
    for t=1:T
        Miss = isnan(X(t,:));
        Obs = ~Miss;
        c=zeros(N,N);
        y = X(t,:)';
        if ~isempty(X(Miss))
            y(Miss) = M(Miss)+S(Miss,Obs)*inv(S(Obs,Obs))*(y(Obs)-M(Obs));
            c(Miss,Miss) = S(Miss,Miss)-S(Miss,Obs)*inv(S(Obs,Obs))*S(Obs,Miss);
        end
	    Y(t,:) = y';
        C(t,:,:) = c + (y-M)*(y-M)';
    end

    % Step 2: update
    M_new = mean(Y)';
    S_new = squeeze(mean(C,1));
    
    D4=[(M_new - M).^4;
                diag((S_new - S).^2) ];
    Distance= mean(D4.^(1/4));
    Convergence = (Distance < Tolerance);
    
    M=M_new;
    S=S_new;
end

E_EM=M;
S_EM=S;