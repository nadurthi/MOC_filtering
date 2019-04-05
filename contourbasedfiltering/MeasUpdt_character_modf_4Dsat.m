function [Xpost_resample,probsXpost_resample,pdfXpostnorm]=MeasUpdt_character_modf_4Dsat(X,probs,pdfnormprior,Nm,Tk,z,Xtruth,model,Xmc,Npts)
% the filter is implemented always using discrete - discrete models
[N,dim]=size(X);

logprobs = log(probs);
% normpdfXmodf = normpdfX;
% normpdfXmodf.func = @(xtrue)eval_truepdf_from_normpdf(xtrue,normpdfX);

% pz1 = integrate_func_exppdf_givenX(@(x)measmodeleval(z,x,model),normpdfXmodf,X,[],[],'GMM_MC');

% Xn = pdfnorm.transForms.trueX2normX(X) ;
% pn = pdfnorm.transForms.trueprob2normprob(probs) ;
% keyboard

% pz2 = integratorFuncTrueX_usingpdfnorm(pdfnormprior,@(x)mvnpdf(repmat(z(:)',size(x,1),1),model.h(x),model.R),'RegTreeBoxIntegrator');

% pz2=1;
% logpz = log(pz2);
pz2 = integratorFuncTrueX_usingpdfnorm(pdfnormprior,@(x)mvnpdf(repmat(z(:)',size(x,1),1),model.hvec(x),model.R),'RegTreeBoxIntegrator');
logpz = log(pz2);



logprobsXpost = zeros(size(probs));
for i=1:size(X,1)
    logprobsXpost(i) = log(1/sqrt(det(2*pi*model.R)))-0.5*(z(:)-model.h(X(i,:)'))'*inv(model.R)*(z(:)-model.h(X(i,:)'))+logprobs(i)-logpz;
end

probsXpost = exp(logprobsXpost);

%% Re-Samplers


[Xpost,postprobs] = SimpleGHpostSampler(X,probs,probsXpost,pdfnormprior,model,z,11,0.5);


% keyboard
%%

figure(103)
plot3(X(:,1),X(:,2),probs/sum(probs),'ro',X(:,1),X(:,2),probsXpost/sum(probsXpost),'gs',Xpost(:,1),Xpost(:,2),postprobs/sum(postprobs),'b+')

% keyboard
%% interpolate pdf probabiliteis at the resamples Xpost points


[mX,PX] = MeanCov(Xpost,postprobs/sum(postprobs));
pdfXpostnorm = get_interp_pdf_0I_4Dsat(Xpost,postprobs,mX,PX,Nm,3,Tk,Xmc,Xtruth);

y=pdfXpostnorm.transForms.trueX2normX(Xpost);
py=pdfXpostnorm.func(y);
probsXpost2=pdfXpostnorm.transForms.normprob2trueprob(py);

Xpost_resample = Xpost;
probsXpost_resample = probsXpost2;

%% Re-sample/ regenerate points
% [~,I] = sort(probsXpost,'descend');
% probsXpost=probsXpost(I);
% X=X(I,:);

% [a,ind]=max(probsXpost);
% mX=X(1,:)';
% PX=zeros(dim);
% wts = probsXpost/sum(probsXpost);
% for i=1:N
%    PX = PX + wts(i)*(X(i,:)'-mX)*(X(i,:)'-mX)';
% end
% [mX,PX] = MeanCov(X(1:50,:),probsXpost(1:50)/sum(probsXpost(1:50)));
% for jj=linspace(3,0.5,10)
%     [Xpost_resample,wsample] = model.pointGenerator(mX,(jj)^2*PX);
%     y=pdfXpostnorm.transForms.trueX2normX(Xpost_resample);
%     py=pdfXpostnorm.func(y);
%     if max(py)/min(py)<500
%         break
%     end
% %     Z=zeros(size(Xpost_resample,1),model.hn);
% %     for i=1:size(Xpost_resample,1)
% %        Z(i,:)=model.h(Xpost_resample(i,:)); 
% %     end
% %     [mz,Pz] = MeanCov(Z,wsample);
% %     Pz=Pz;
% %     pzk = mvnpdf(z(:)',mz(:)',Pz);
% %     pzm = mvnpdf(mz(:)',mz(:)',Pz);
% %     if pzk>pzm/10
% %         break
% %     end
%     jj
% end
% keyboard



% [Xpost_resample,wsample] = model.pointGenerator(mX,(0.5)^2*PX);
% Xpost_resample=X;

% y=pdfXpostnorm.transForms.trueX2normX(Xpost_resample);
% py=pdfXpostnorm.func(y);
% probsXpost_resample=pdfXpostnorm.transForms.normprob2trueprob(py);



