function [Xpost,postprobs] = SimpleMCpostSampler(Xprior,priorprobs,probsXpost,priorpdfnorm,model,z,N1d,factor)
% Just compute posterior mX,PX and get the samples from that.
% priorprobs,probsXpost are both at the points Xprior
% N1d is the number of points in 1D

[N,dim]=size(Xprior);

%% Get New prior points
[mX,PX] = MeanCov(Xprior,probsXpost/sum(probsXpost));
[XpriorNew,~] = model.pointGenerator(mX,(factor)^2*PX);

%% get new prior point's probnabilities
y=priorpdfnorm.transForms.trueX2normX(XpriorNew);
py=priorpdfnorm.func(y);
priorprobsNew=priorpdfnorm.transForms.normprob2trueprob(py);
logprobs = log(priorprobsNew);

%% compute bayes constants using prior pdfnorm
% keyboard
pz2 = integratorFuncTrueX_usingpdfnorm(priorpdfnorm,@(x)mvnpdf(repmat(z(:)',size(x,1),1),model.hvec(x),model.R),'RegTreeBoxIntegrator');
logpz = log(pz2);

% tic
% priorpdfnorm.func(Xprior(1:200,:))
% toc

%% update the new prior points
logprobsXpost = zeros(size(logprobs));
for i=1:size(XpriorNew,1)
    logprobsXpost(i) = log(1/sqrt(det(2*pi*model.R)))-0.5*(z(:)-model.h(XpriorNew(i,:)'))'*inv(model.R)*(z(:)-model.h(XpriorNew(i,:)'))+logprobs(i)-logpz;
end

Xpost = XpriorNew;
postprobs = exp(logprobsXpost);