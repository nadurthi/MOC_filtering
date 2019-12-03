function [pnorm,logpnorm,xk,xknorm] = getbacknormprobs(model,Tk1,dtkk1,normpdf_atk,normpdfexp_atk,xnorm,transForms_atk1,transForms_atk)


xk1=transForms_atk1.normX2trueX(xnorm) ;

xk = model.fback(dtkk1,Tk1,xk1);
xknorm = transForms_atk.trueX2normX(xk);

%% eval at norm pdf at prev time k
pknorm = normpdf_atk(xknorm);
logpknorm = normpdfexp_atk(xknorm);

% transForms.trueprob2normprob_constant
% transForms.normprob2trueprob_constant
%% now trqansform from norm to true space at time k 
pk = transForms_atk.normprob2trueprob(pknorm);
logpk = log(transForms_atk.normprob2trueprob_constant)+logpknorm;

%% prop from k to k1
pk1 = pk; % as sat dynamics do not diverge
logpk1 = logpk;

%% transform from true space to norm space at time k1
pnorm=transForms_atk1.trueprob2normprob(pk1);
logpnorm = log(transForms_atk1.trueprob2normprob_constant)+logpk1;


end