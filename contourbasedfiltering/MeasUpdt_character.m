function [Xpost_resample,probsXpost_resample]=MeasUpdt_character(normpdfX,X,probs,mquad,Pquad,Nm,Tk,z,model)
% the filter is implemented always using discrete - discrete models
logprobs = log(probs);
normpdfXmodf = normpdfX;
normpdfXmodf.func = @(xtrue)eval_exppdf_from_normpdf(xtrue,normpdfX);

pz1 = integrate_func_exppdf_givenX(@(x)gausspdfmodeleval(z,x,model.R),normpdfXmodf,X,[],[],'GMM_MC');

Xn=zeros(size(X));
for i=1:size(X,1)
   Xn(i,:)=normpdfX.trueX2normX(X(i,:)) ;
end
pz2 = integrate_func_exppdf_givenX(@(x)gausspdfmodeleval(z,normpdfX.normX2trueX(x),model.R),normpdfX,Xn,[],[],'GMM_MC');

logpz = log(pz1);

logprobsXpost = zeros(size(probs));
for i=1:size(X,1)
    logprobsXpost(i) = log(1/sqrt(det(2*pi*model.R)))-0.5*(z(:)-model.h(X(i,:)'))'*inv(model.R)*(z(:)-model.h(X(i,:)'))+logprobs(i)-logpz;
end
probsXpost = exp(logprobsXpost);


%% Estimate normalizing constant

pdfXpostnorm = get_interp_pdf_0I(X,probsXpost,mquad,Pquad,Nm,Tk,[]);
y=pdfXpostnorm.trueX2normY(X);
py=pdfXpostnorm.func(y);
probsXpost2=pdfXpostnorm.normprob2trueprob(py);
%% Re-sample/ regenerate points

[Xpost_resample,~] = GH_points(mquad,0.5^2*Pquad,5);

y=pdfXpostnorm.trueX2normY(Xpost_resample);
py=pdfXpostnorm.func(y);
probsXpost_resample=pdfXpostnorm.normprob2trueprob(py);



