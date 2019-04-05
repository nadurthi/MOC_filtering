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

figure(111)
plot3(Xprior_norm(:,1),Xprior_norm(:,2),priorprobs,'ro',Xprior_norm(:,1),Xprior_norm(:,2),probsXpost,'b+')


% GMMHull.GMM=[];
% GMMHull.FitGMM_dummy(Xprior_norm,15,0.001)
GMMHull.resetGMM_Hull_full()
GMMHull.SetGMM_Hull(12,0.01);
GMMHull.GMM=GMMHull.GMMhull;
figure(114)
GMMHull.plotGMMpoints([1,2],'ro')
figure(115)
GMMHull.plotGMMpoints([3,4],'ro')



% GMMHull.optimGMMhullwts_relative2probs(Xprior_norm,priorprobs_norm);
GMMHull.optimGMMhullwts_reoptimize(Xprior_norm,postprobs_norm);
GMMHull.GMM=GMMHull.GMMhull;


figure(113)
states=[1,2];
GMMHull.plotGMMSurf(states,'g')
hold on
plot3(Xprior_norm(:,states(1)),Xprior_norm(:,states(2)),priorprobs_norm/max(priorprobs_norm),'ro')
hold off
figure(118)
states=[3,4];
GMMHull.plotGMMSurf(states,'g')
hold on
plot3(Xprior_norm(:,states(1)),Xprior_norm(:,states(2)),priorprobs_norm/max(priorprobs_norm),'ro')
hold off




% GMMHull.optimGMMwts_relative2probs_and_setmXPx(Xprior_norm,postprobs_norm);
if dim==2
    XpriorNew_norm=GMMHull.gen_quad_GMM_GH_2D(11);
elseif dim==4
    XpriorNew_norm=GMMHull.gen_quad_GMM_GH_4D(11);
end
XpriorNew = priorpdfnorm.transForms.normX2trueX(XpriorNew_norm);

XpriorNew=[XpriorNew;Xprior(probsXpost>max(probsXpost)/10,:)];





% XpriorNew = Xprior;
%% get new prior point's probnabilities
XpriorNew_norm=priorpdfnorm.transForms.trueX2normX(XpriorNew);
priorprobsNew_norm=priorpdfnorm.func(XpriorNew_norm);
priorprobsNew=priorpdfnorm.transForms.normprob2trueprob(priorprobsNew_norm);



% keyboard
% 
% XpriorNew=XpriorNew(priorprobsNew>(max(priorprobsNew)/10),:);
% priorprobsNew=priorprobsNew(priorprobsNew>(max(priorprobsNew)/10));


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

XpriorNew_beforereduction = XpriorNew;

if length(postprobs)>7000
    frac=7000/length(postprobs);
    qq=quantile(postprobs,1-frac);
    ind = postprobs>qq;
    postprobs = postprobs(ind);
    priorprobsNew_reduced = priorprobsNew(ind);
    Xpost = Xpost(ind,:);
    XpriorNew = XpriorNew(ind,:);
    
    
    
end



%%
XpriorNew_norm=priorpdfnorm.transForms.trueX2normX(XpriorNew);
priorprobsNew_norm=priorpdfnorm.func(XpriorNew_norm);

XpriorNew_beforereduction_norm=priorpdfnorm.transForms.trueX2normX(XpriorNew_beforereduction);


postprobsNew_norm=priorpdfnorm.transForms.trueprob2normprob(postprobs);


figure(119)
plot3(Xprior_norm(:,1),Xprior_norm(:,2),priorprobs_norm/sum(priorprobs_norm),'ro',Xprior_norm(:,1),Xprior_norm(:,2),postprobs_norm/sum(postprobs_norm),'b+',XpriorNew_norm(:,1),XpriorNew_norm(:,2),priorprobsNew_norm/sum(priorprobsNew_norm),'gs',XpriorNew_norm(:,1),XpriorNew_norm(:,2),postprobsNew_norm/sum(postprobsNew_norm),'k^')
legend('prior','priorpt-postprob','newpts-priorprobs','newpts-postprobs')



% keyboard

figure(112)
plot3(Xprior(:,1),Xprior(:,2),priorprobs/sum(priorprobs),'ro',Xprior(:,1),Xprior(:,2),probsXpost/sum(probsXpost),'b+',XpriorNew(:,1),XpriorNew(:,2),priorprobsNew_reduced/sum(priorprobsNew_reduced),'gs',XpriorNew(:,1),XpriorNew(:,2),postprobs/sum(postprobs),'k^')
legend('prior','priorpt-postprob','newpts-priorprobs','newpts-newpostprobs')

figure(113)
GMMHull.plotGMMSurf([1,2],'g')
hold on
y=priorpdfnorm.transForms.trueX2normX(Xpost);
plot(y(:,1),y(:,2),'k*')
hold off
Fpost=GMMHull.GMMhull.w

% keyboard


