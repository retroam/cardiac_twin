function error = objectiveFcn(KL,KR,gL,aA,gA,paramsEst,paramsFixed,simOptions,tspan,y0,data)
% objective function for performing least squares minimization
% This error function is weighted by the std dev of the experimental
% measurements, which affects confidence intervals and relative weighting 
% of multi-objective function estimates.

params = paramsFixed;

gamma_L = gL; %1
alpha_A = aA; %1.5574
gamma_A = gA; %1 
alpha_L=paramsEst;
params(8)  = (alpha_L*KL); 
params(12) = (alpha_L*KR);
params(16) = (alpha_L*gamma_L*KL);
params(20)  = (alpha_L*KL);
params(24) = (alpha_A*KR);
params(26) = (alpha_A*gamma_A*KL);

[tSim,ySim] = ode15s(@daeODE,tspan,y0,simOptions,params);    % run the model for a given sample of paramsEst
cAMP = ySim(end,9);
error = cAMP - data; 
