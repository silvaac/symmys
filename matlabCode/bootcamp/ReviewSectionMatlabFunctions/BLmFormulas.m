function [BLMu,BLSigma]=BLmFormulas(Mu,Sigma,P,v,Omega)

BLMu = Mu + Sigma*P'*inv(P*Sigma*P'+Omega)*(v-P*Mu);
BLSigma =  Sigma -  Sigma*P'*inv(P*Sigma*P'+Omega)*P*Sigma;
