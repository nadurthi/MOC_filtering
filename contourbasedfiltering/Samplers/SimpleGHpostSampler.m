function [Xpost,postprobs] = SimpleGHpostSampler(Xprior,priorprobs,probsXpost,priorpdfnorm,model,z,N1d,factor)
% Just compute posterior mX,PX and get the samples from that.
% priorprobs,probsXpost are both at the points Xprior
% N1d is the number of points in 1D

[N,dim]=size(Xprior);

%% Get New prior points
[mXprior,PXprior] = MeanCov(Xprior,priorprobs/sum(priorprobs));

[mXpost,PX] = MeanCov(Xprior,probsXpost/sum(probsXpost));

lambda =0.9;
mX = lambda*mXpost+(1-lambda)*mXprior;


[XpriorNew,~] = model.pointGenerator(mX,(factor)^2*PX);

% keyboard

GMMHull=priorpdfnorm.GMMHull;
priorprobs_norm = priorpdfnorm.transForms.trueprob2normprob(priorprobs);
postprobs_norm = priorpdfnorm.transForms.trueprob2normprob(probsXpost);
Xprior_norm = priorpdfnorm.transForms.trueX2normX(Xprior);


% GMMHull.GMM=[];
% GMMHull.FitGMM_dummy(Xprior_norm,15,0.001)
GMMHull.resetGMM_Hull_full()
GMMHull.SetGMM_Hull(10,0.01);
GMMHull.GMM=GMMHull.GMMhull;
% GMMHull.optimGMMhullwts_relative2probs(Xprior_norm,postprobs_norm);
GMMHull.optimGMMhullwts_reoptimize(Xprior_norm,postprobs_norm);

% GMMHull.optimGMMwts_relative2probs_and_setmXPx(Xprior_norm,postprobs_norm);

XpriorNew_norm=GMMHull.gen_quad_GMM_GH_2D(11);
XpriorNew = priorpdfnorm.transForms.normX2trueX(XpriorNew_norm);

XpriorNew=[XpriorNew;Xprior(probsXpost>max(probsXpost)/10,:)];



% XpriorNew = Xprior;
%% get new prior point's probnabilities
y=priorpdfnorm.transForms.trueX2normX(XpriorNew);
py=priorpdfnorm.func(y);
priorprobsNew=priorpdfnorm.transForms.normprob2trueprob(py);

XpriorNew=XpriorNew(priorprobsNew>(max(priorprobsNew)/10),:);
priorprobsNew=priorprobsNew(priorprobsNew>(max(priorprobsNew)/10));


logpriorprobsnew = log(priorprobsNew);



%% compute bayes constants using prior pdfnorm
pz2 = integratorFuncTrueX_usingpdfnorm(priorpdfnorm,@(x)mvnpdf(repmat(z(:)',size(x,1),1),model.hvec(x),model.R),'RegTreeBoxIntegrator');
logpz = log(pz2);

%% update the new prior points
logprobsXpost = zeros(size(logpriorprobsnew));
for i=1:size(XpriorNew,1)
    logprobsXpost(i) = log(1/sqrt(det(2*pi*model.R)))-0.5*(z(:)-model.h(XpriorNew(i,:)'))'*inv(model.R)*(z(:)-model.h(XpriorNew(i,:)'))+logpriorprobsnew(i)-logpz;
end

Xpost = XpriorNew;
postprobs = exp(logprobsXpost);

%%
figure(111)
plot3(Xprior_norm(:,1),Xprior_norm(:,2),priorprobs,'ro',Xprior_norm(:,1),Xprior_norm(:,2),probsXpost,'b+')

figure(112)
plot3(Xprior(:,1),Xprior(:,2),priorprobs/sum(priorprobs),'ro',Xprior(:,1),Xprior(:,2),probsXpost/sum(probsXpost),'b+',XpriorNew(:,1),XpriorNew(:,2),priorprobsNew/sum(priorprobsNew),'gs',XpriorNew(:,1),XpriorNew(:,2),postprobs/sum(postprobs),'k^')
legend('prior','priorpt-postprob','newpts-priorprobs','newpts-newpostprobs')

figure(113)
GMMHull.plotGMMSurf([1,2],'g')
hold on
y=priorpdfnorm.transForms.trueX2normX(Xpost);
plot(y(:,1),y(:,2),'k*')
hold off
Fpost=GMMHull.GMMhull.w

% keyboard


