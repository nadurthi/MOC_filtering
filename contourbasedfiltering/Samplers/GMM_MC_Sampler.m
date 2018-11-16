function GMM_MC_Sampler()
gmm = MyGmm2MatlabGMM(pdfnormprior.GMMHull.GMMhull);
XmcHull = random(gmm,10000);
indbnd = pdfnormprior.GMMHull.IsInsideHull(XmcHull,1.0);
Y=zeros(10000,1);
Y(indbnd)=100;
% XmcHull = XmcHull(indbnd,:);
Nmc = size(XmcHull,1);

figure(104)
Xn = pdfnormprior.transForms.trueX2normX(X);
Xnpost = pdfnormprior.transForms.trueX2normX(Xpost);
% pn = pdfnormprior.transForms.normprob2trueprob(probs);
plot3(Xn(:,1),Xn(:,2),probs/sum(probs),'ro',Xn(:,1),Xn(:,2),probsXpost/sum(probsXpost),'b+')
hold on
plot3(XmcHull(:,1),XmcHull(:,2),-1*ones(Nmc,1),'.')
plot(Xnpost(:,1),Xnpost(:,2),'gs')
pdfnormprior.GMMHull.plotGMMpointsHUll([1,2],[],1.2,'ro')
hold off