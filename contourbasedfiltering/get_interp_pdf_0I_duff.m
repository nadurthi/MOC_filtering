function pdfnorm = get_interp_pdf_0I_duff(X,probs,mquad,Pquad,Nm,Ngcomp,Tk,Xmctest,Xtruth)
%%
[N,dim] =size(X);

% probs = PROBS;
% X=XX;
% mquad = MMQ;
% Pquad = PPQUAD;
% 
% sqP = sqrtm(Pquad);
% diag_sqP = diag(sqP);
% bndslb = mquad(:)'-3*diag_sqP(:)';
% bndsub = mquad(:)'+3*diag_sqP(:)';
% 
% [X,probs]=filterX_inBox(X,probs,bndslb,bndsub);
% [mquad,Pquad]=MeanCov(X,probs/sum(probs));

dsX = DataSet(X,probs,'TrueState');


mquad=mquad(:);

dsX.AddMeanCov_to_OI_Trasform(mquad,Pquad);
dsX.AddHyperCubeTrasform(-1*ones(dim,1),1*ones(dim,1));

% plot3(dsX.X(:,1),dsX.X(:,2),dsX.p,'ro')

% keyboard

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

% GMMfitter = GMMFitDataSet(dsX.X,dsX.p);
% % GMM = GMMfitter.FitGMM_1comp();
% GMM = GMMfitter.FitGMM_kmeans_optimwt(5);
% 
% % figure(36)
% % GMMfitter.plotGMMpointsHUll([1,2],2,'ro')
% 
% figure(34)
% GMMfitter.plotGMMpoints([1,2],'ro')
% title('GMM fit points')
% hold off
% 
% figure(35)
% GMMfitter.plotGMMSurf([1,2],'ro')
% title('GMM fit surf')
% hold off

% keyboard
%% fitting poly to log of probas
% keyboard

Pf=Basis_polyND(dim,Nm);



%% reguralization points

Xineq=[];
[Xbnd1,~] = GLgn_pts(-1.5*ones(1,dim),1.5*ones(1,dim),9);

Xbnd2=[];
% [Xbnd2,~] = GLgn_pts(-2*ones(1,dim),2*ones(1,dim),9);


XtestoutsideHull = mvurnd(-1.7*ones(dim,1),1.7*ones(dim,1),3000);

% Xbnd2=2.5*gen_uniform_grid(3,dim);
% Xbnd2=5*gen_uniform_grid(5,dim);
Xineq = [Xbnd1;Xbnd2];

% keyboard
GMMHull = GMMFitDataSet(dsX.X,dsX.p);
GMMHull.SetGMM_Hull(10);
% GMMHull.optimGMMhullwts();
% GMMHull.GMM=GMMHull.GMMhull;


% figure
% GMMHull.plotGMMSurf([1,2],'g')

indbnd = GMMHull.IsInsideHull(Xineq,1.5);
Xineq = Xineq(~indbnd,:);

indbnd = GMMHull.IsInsideHull(XtestoutsideHull,1.5);
XtestoutsideHull = XtestoutsideHull(~indbnd,:);


figure(36)
GMMHull.plotGMMpointsHUll([1,2],Xineq,1.0,'ro')
hold on
plot3(Xineq(:,1),Xineq(:,2),-ones(size(Xineq,1),1),'b+')
hold off
title('GMM HUll')

Nineq = size(Xineq,1);
%% full pdf without gaussian
% keyboard

LB = -2*ones(dim,1);
UB = 2*ones(dim,1);

LBtest=-1.5*ones(dim,1);
UBtest=1.5*ones(dim,1);

% XtestingMC = mvurnd(LB,UB,20000);
    
    RIF=RegInterpFitters('DecisionTreeAdaptiveOutRegion');
    RIF.fit(dsX.X,dsX.p,Xineq,[],[],GMMHull)
    
    RigTree = RIF;
    
% while 1

    RIF=RegInterpFitters('PolyFit');%'KnnMean''PolyFit'
    RIF.fit(dsX.X,dsX.p,Xineq,XtestoutsideHull,Pf,GMMHull);

%     RIF.plot(dsX.X,dsX.p,Xineq,[1,2],-1.5*ones(1,dim),1.5*ones(1,dim))

%     ExpPolyfitter = PolyFit(dsX.X,dsX.p);
%     mxentpoly_norm =ExpPolyfitter.fitExpPoly_A_Atop_Aineq(Pf,Xineq,XtestingMC);
%     mxentpoly_norm =ExpPolyfitter.fitExpPoly_A_Atop_AdaptiveC(Pf,XtestingMC);
    
    SS = mvurnd(LB,UB,10000);
    [pSS,plogSS]=RIF.evalfit(SS);
    
%     pSS = evaluate_PolyforLargetSetX(mxentpoly_norm,SS);
    ind = pSS>max(log(dsX.p));
    
    figure(37)
    states=[1,2];
    RIF.plot(dsX.X,dsX.p,SS,states,LB,UB);
    hold off
    
%     figure(40)
%     ExpPolyfitter.PlotExpPolyFits_points([1,2],-1*ones(dim,1),ones(dim,1))
%     hold off
%     
    
%     keyboard
    
%     if sum(ind)==0
%         disp('All probs are in the bounds')
%         break
%     end
    % Xbnd2=2.5*gen_uniform_grid(3,dim);
%     Xineq = [Xineq;SS(ind,:)];
    
    

% end
% fitstats_norm = TestPolyFits(mxentpoly_norm,dsX.X,log(dsX.p),LBtest,UBtest);

if isfield(RIF.method_params,'mxentpoly')
    pdfnorm.dim =dim;
    mxentpoly_norm = RIF.method_params.mxentpoly;
    pdfnorm.mxentpoly_norm = mxentpoly_norm;
    pdfnorm.func=@(x)exp(evaluate_polyND(mxentpoly_norm,x));
    pdfnorm.transForms = dsX.GetTrasnformers();
    
    pdfnorm.info = 'true-0I-hypercube-11';
    pdfnorm.pdftype = 'ExpPdf';
    
    pdfnorm.GMMHull = GMMHull;
    pdfnorm.LB = LB;
    pdfnorm.UB = UB;
    pdfnorm.RigTree = RigTree;
    
    normconst = integratorFuncTrueX_usingpdfnorm(pdfnorm,@(x)constantfunc(x,1),'RegTreeBoxIntegrator');
%     normconst=1;
    pdfnorm.func=@(x)(1/normconst)*exp(evaluate_polyND(mxentpoly_norm,x));
    
%     pdfnorm.polyeval=@(x)evaluate_polyND(mxentpoly_norm,x);
%     pdfnorm.poly=mxentpoly_norm;
    
    
else
    mxentpoly_norm=[];
    
    pdfnorm.dim =dim;
    pdfnorm.func=@(x)RIF.evalfit(x);
    pdfnorm.transForms = dsX.GetTrasnformers();
    
    pdfnorm.info = 'true-0I-hypercube-11';
    pdfnorm.pdftype = 'HybridPdf';
    pdfnorm.GMMHull = GMMHull;
    pdfnorm.LB = LB;
    pdfnorm.UB = UB;
    pdfnorm.RigTree = RigTree;
    
    normconst = integratorFuncTrueX_usingpdfnorm(pdfnorm,@(x)constantfunc(x,1),'RegTreeBoxIntegrator');
    pdfnorm.func=@(x)(1/normconst)*RIF.evalfit(x);
    
    
    
end

% keyboard


% pdfnorm = normalize_exp_pdf(pdfnorm,dsX.X,dsX.p,mquad,Pquad,[],'GMM_MC');

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

