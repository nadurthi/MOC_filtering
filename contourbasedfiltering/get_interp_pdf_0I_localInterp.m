function pdfnorm = get_interp_pdf_0I_localInterp(X,probs,mquad,Pquad,Nm,Ngcomp,Tk,Xmctest,Xtruth)
%%
dsX = DataSet(X,probs,'TrueState');
[N,dim] =size(X);

mquad=mquad(:);

dsX.AddMeanCov_to_OI_Trasform(mquad,Pquad);
dsX.AddHyperCubeTrasform(-1*ones(dim,1),1*ones(dim,1));

indd = dsX.p<1e-70;
dsX.p(indd) = 1e-70;

dsX.SortByProb('descend');

logpn = dsX.getlogProb();

if isempty(Xmctest)==0
    [Xnmctest,~] = dsX.ApplyAffineTransform_Original2Final(Xmctest, zeros(size(Xmctest,1),1) );
end
if isempty(Xtruth)==0
    [Xntruth,~] =  dsX.ApplyAffineTransform_Original2Final(Xtruth, zeros(size(Xtruth,1),1) );
end


%% plottinmg

figure(33)
dsX.PlotPointsProbs3D([1,2],'ro');
hold on 
if isempty(Xmctest)==0
plot(Xnmctest(:,1),Xnmctest(:,2),'bs')
end
if isempty(Xtruth)==0
plot(Xntruth(:,1),Xntruth(:,2),'k*')
end
title(['true points and MC: time step = ',num2str(Tk)])
hold off

%% First fit Ngcomp gaussian components to the points

GMMfitter = GMMFitDataSet(dsX.X,dsX.p);
GMM = GMMfitter.FitGMM_kmeans_optimwt(5);

figure(34)
GMMfitter.plotGMMpoints([1,2],'ro')
title('GMM fit points')
hold off

figure(35)
GMMfitter.plotGMMSurf([1,2],'ro')
title('GMM fit surf')
hold off


%% fitting poly to log of probas
% keyboard

% Pf=Basis_polyND(dim,Nm);



%% reguralization points

Xineq=[];
% [Xbnd1,~] = GLgn_pts(-1.5*ones(1,dim),1.5*ones(1,dim),4);

SS = mvurnd(-1.5*ones(dim,1),1.5*ones(dim,1),50000);

% Xbnd2=2.5*gen_uniform_grid(3,dim);
% Xbnd2=5*gen_uniform_grid(5,dim);
Xineq = [SS];

GMMfitter.SetGMM_Hull(6);
indbnd = GMMfitter.IsInsideHull(Xineq,1);
Xineq = Xineq(~indbnd,:);

figure(36)
GMMfitter.plotGMMpointsHUll([1,2],1,'ro')
hold on
plot3(Xineq(:,1),Xineq(:,2),-ones(size(Xineq,1),1),'b+')
hold off
title('GMM HUll')

Nineq = size(Xineq,1);
%% full pdf without gaussian
% keyboard

LB=-1.5*ones(dim,1);
UB=1.5*ones(dim,1);

LBtest=-1.5*ones(dim,1);
UBtest=1.5*ones(dim,1);



% fitstats_norm = TestPolyFits(mxentpoly_norm,dsX.X,log(dsX.p),LBtest,UBtest);

pdfnorm.func=@(x)estimateprobknn(x,[dsX.X;Xineq],[dsX.p;1e-3*min(dsX.p)*ones(Nineq,1)],3);
pdfnorm.polyeval=NaN;
pdfnorm.poly=NaN;
pdfnorm.info = 'knn';
pdfnorm.pdftype = 'HybridPdf';

% keyboard

pdfnorm = normalize_exp_pdf(pdfnorm,dsX.X,dsX.p,mquad,Pquad,GMM,'GMM_MC');
pdfnorm.transForms = dsX.GetTrasnformers();
% [Atranf2norm,mtransf2norm]=obj.GetAffineTransform_Original2Final();
% [Atranf2true,mtransf2true]=obj.GetAffineTransform_Final2Original();
% transForms.trueX2normX = @(x)affineTransform(x,Atranf2norm,mtransf2norm);
% transForms.normX2trueX = @(xn)affineTransform(xn,Atranf2true,mtransf2true);
% transForms.normprob2trueprob = @(p)p/det(Atranf2true);
% transForms.trueprob2normprob = @(p)p/det(Atranf2norm);


%%
% cm=1;
% pdiff=(dsX.p)./GaussSumMix(dsX.X,GMM);
% if any(pdiff<=0)
%     keyboard
% end
% 
% % mndiff = min(diff) - cm;
% % mndiff = 1*mndiff;
% % pdiff = diff - mndiff;
% 
% % Pf=Basis_polyND(dim,6);
% ExpPolyfitter = PolyFit(dsX.X,pdiff);
% mxentpoly_normGM = ExpPolyfitter.fitExpPoly_A_Atop_Aineq(Pf,Xineq);
% figure(38)
% ExpPolyfitter.PlotExpPolyFits([1,2],LB,UB)
% title('gaussian base: diff fit')
% hold off
% 
% fitstats_normGM = TestPolyFits(mxentpoly_normGM,dsX.X,log(pdiff),LBtest,UBtest);
% 
% pdfnormGM1.func =@(x)GaussSumMix(x,GMM).*exp(evaluate_polyND(mxentpoly_normGM,x));
% pdfnormGM1.info = 'SingleGMMComp-true-0I-hypercube-11';
% pdfnormGM1.pdftype = 'HybridPdf';
% 
% pdfnormGM1 = normalize_exp_pdf(pdfnormGM1,dsX.X,dsX.p,mquad,Pquad,'GMM_MC');
% pdfnormGM1.transForms = dsX.GetTrasnformers();
% 
% keyboard

%% Fit statistics

% ErrApproaches = [100*(pdfnormGM1.func(dsX.X)-dsX.p)./dsX.p,100*(pdfnorm.func(dsX.X)-dsX.p)./dsX.p];
ErrApproaches = [100*(log(pdfnorm.func(dsX.X))-log(dsX.p))./log(dsX.p)];
mmbin = min(min(ErrApproaches,[],1));
mxbin = max(max(ErrApproaches,[],1));

mmbin=max(mmbin,-300);
mxbin=min(mxbin,300);

figure(39)
histogram(ErrApproaches(:,1),100,'facecolor','r','facealpha',0.5) 
% hold on
% histogram(ErrApproaches(:,2),linspace(mmbin,mxbin,100),'facecolor','b','facealpha',0.5) 
legend('Full')
% hold off
title(['Error st for k = ',num2str(Tk)])

% keyboard

