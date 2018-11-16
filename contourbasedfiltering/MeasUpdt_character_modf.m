function [Xpost_resample,probsXpost_resample,pdfXpostnorm]=MeasUpdt_character_modf(X,probs,Nm,Tk,z,Xtruth,model,Xmc,Npts)
% the filter is implemented always using discrete - discrete models


logprobs = log(probs);
% normpdfXmodf = normpdfX;
% normpdfXmodf.func = @(xtrue)eval_truepdf_from_normpdf(xtrue,normpdfX);

% pz1 = integrate_func_exppdf_givenX(@(x)measmodeleval(z,x,model),normpdfXmodf,X,[],[],'GMM_MC');

% Xn = pdfnorm.transForms.trueX2normX(X) ;
% pn = pdfnorm.transForms.trueprob2normprob(probs) ;
% keyboard
% pz2 = integrate_func_exppdf_givenX(@(x)measmodeleval(z,normpdfX.normX2trueX(x),model),normpdfX,Xn,[],[],'GMM_MC');

pz2=1;
logpz = log(pz2);


logprobsXpost = zeros(size(probs));
for i=1:size(X,1)
    logprobsXpost(i) = log(1/sqrt(det(2*pi*model.R)))-0.5*(z(:)-model.h(X(i,:)'))'*inv(model.R)*(z(:)-model.h(X(i,:)'))+logprobs(i)-logpz;
end
probsXpost = exp(logprobsXpost);

% figure(8)
% plot3(X(:,1),X(:,1),probs,'ro',X(:,1),X(:,1),probsXpost,'b+')

%% Estimate normalizing constant

% keyboard
figure
plot3(X(:,1),X(:,2),probs,'ro')
hold on
plot3(X(:,1),X(:,2),probsXpost,'b+')

[mX,PX] = MeanCov(X,probsXpost/sum(probsXpost));
% fullnormpdf=get_interp_pdf_0I_boostmixGaussian(X,probs,mX,PX,4,3,k,Xmctest,Xtruth(k,:)); %Xtruth(k,:)
% pdfXpostnorm = get_interp_pdf_0I_boostmixGaussian(X,probsXpost,mX,PX,Nm,3,Tk,Xmc,Xtruth);
pdfXpostnorm = get_interp_pdf_0I_duff(X,probsXpost,mX,PX,Nm,3,Tk,Xmc,Xtruth);
% y=pdfXpostnorm.trueX2normX(X);
% py=pdfXpostnorm.func(y);
% probsXpost2=pdfXpostnorm.normprob2trueprob(py);

%% Re-sample/ regenerate points
[Xpost_resample,~] = model.pointGenerator(mX,0.5^2*PX);
% [Xpost_resample,~] = GH_points(mX,0.5^2*PX,Npts);

y=pdfXpostnorm.transForms.trueX2normX(Xpost_resample);
py=pdfXpostnorm.func(y);
probsXpost_resample=pdfXpostnorm.transForms.normprob2trueprob(py);



