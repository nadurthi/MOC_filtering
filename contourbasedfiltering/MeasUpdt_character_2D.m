function [Xpost_resample,probsXpost_resample,Xquad_post,wquad_post,pdfXpostnorm]=MeasUpdt_character_2D(normpdfX,X,probs,Xquad,wquad,Nm,Tk,z,model,Xmc)
% the filter is implemented always using discrete - discrete models


logprobs = log(probs);
normpdfXmodf = normpdfX;
normpdfXmodf.func = @(xtrue)eval_exppdf_from_normpdf(xtrue,normpdfX);

% pz1 = integrate_func_exppdf_givenX(@(x)measmodeleval(z,x,model),normpdfXmodf,X,[],[],'GMM_MC');

Xn=normpdfX.trueX2normX(X) ;

% keyboard
pz2 = integrate_func_exppdf_givenX(@(x)measmodeleval(z,normpdfX.normX2trueX(x),model),normpdfX,Xn,[],[],'GMM_MC');

logpz = log(pz2);


logprobsXpost = zeros(size(probs));
for i=1:size(X,1)
    logprobsXpost(i) = log(1/sqrt(det(2*pi*model.R)))-0.5*(z(:)-model.h(X(i,:)'))'*inv(model.R)*(z(:)-model.h(X(i,:)'))+logprobs(i)-logpz;
end
probsXpost = exp(logprobsXpost);

% figure(8)
% plot3(X(:,1),X(:,1),probs,'ro',X(:,1),X(:,1),probsXpost,'b+')
%% Run ukf
% [Xq,wq]= model.quadfunc(mquad,Pquad);
Xq=Xquad;
wq=wquad;
[mquad,Pquad] = MeanCov(Xq,wq);

Zq=zeros(size(Xq,1),model.hn);
for i=1:size(Xq,1)
   Zq(i,:) = model.h(Xq(i,:)); 
end
[mz,Pz] = MeanCov(Zq,wq);
Pz = Pz+model.R;
Pcc=CrossCov(Xq,mquad,Zq,mz,wq);
K=Pcc/Pz;
mquad_post=mquad+K*(z(:)-mz(:));
Pquad_post=Pquad-K*Pz*K';
[Xquad_post,wquad_post]= model.quadfunc(mquad_post,Pquad_post);
%% Estimate normalizing constant

[mquad_post,Pquad_post] = MeanCov(X,probsXpost/sum(probsXpost));

pdfXpostnorm = get_interp_pdf_0I(X,probsXpost,mquad_post,Pquad_post,Nm,Tk,Xmc);
y=pdfXpostnorm.trueX2normX(X);
py=pdfXpostnorm.func(y);
probsXpost2=pdfXpostnorm.normprob2trueprob(py);

[Xquad_post,wquad_post]= model.quadfunc(mquad_post,Pquad_post);

%% Re-sample/ regenerate points

[Xpost_resample,~] = GH_points(mquad_post,0.5^2*Pquad_post,5);

y=pdfXpostnorm.trueX2normX(Xpost_resample);
py=pdfXpostnorm.func(y);
probsXpost_resample=pdfXpostnorm.normprob2trueprob(py);



