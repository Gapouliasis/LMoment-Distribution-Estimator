function par=gevLMfit(X,bias)
% Calculates the L-Moments estimators of the Generelized Extreme Value Distribution
% with either the biased or unbiased L-Moment estimator 
% Written by George Pouliasis
LM=LMoments(X,bias);
tau3=LM(3)/LM(2);
c=log(2)/log(3)-2/(3+tau3);
k=7.8*c-1.43*c^2;
lamda=k*LM(2)/(gamma(1-k)*(2^k-1));
psi=LM(1)/lamda-(gamma(1-k)-1)/k;
par(1)=k;
par(2)=lamda;
par(3)=psi*lamda;
end