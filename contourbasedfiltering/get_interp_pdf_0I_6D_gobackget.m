function pdfnorm = get_interp_pdf_0I_6D_gobackget(model,Xk,probsk,Nm,Tstepk1,Tk1,dtkk1,priorpdfnormk,Xmctest,Xtruth)
% known prior pdf priorpdfnormk at time k
% compute pdfnorm or pdfnormk1 at time k+1
% Tk1 is time at k+1,
% dtkk1 is deltat from k to k+1
%%
% [N,dim] =size(X);

[mXk1,PXk1]=MeanCov(Xk,probsk/sum(probsk));

mXk1=mXk1(:);
dim = length(mXk1);

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

dsX = DataSet(Xk,probsk,'TrueState');


dsX.AddMeanCov_to_OI_Trasform(mXk1,PXk1);
dsX.AddHyperCubeTrasform(-1*ones(dim,1),1*ones(dim,1));

% plot3(dsX.X(:,1),dsX.X(:,2),dsX.p,'ro')

% keyboard

% indd = dsX.p<1e-100;
% dsX.p(indd) = 1e-100;

dsX.SortByProb('descend');

logpn = dsX.getlogProb();

if isempty(Xmctest)==0
    [Xnmctest,~] = dsX.ApplyAffineTransform_Original2Final(Xmctest, zeros(size(Xmctest,1),1) );
end
if isempty(Xtruth)==0
    [Xntruth,~] =  dsX.ApplyAffineTransform_Original2Final(Xtruth, zeros(size(Xtruth,1),1) );
end


%% plottinmg
for i=0:4
    figure(33+i)
    dsX.PlotPointsProbs3D([i+1,i+2],'ro');
    hold on
    if isempty(Xmctest)==0
        plot(Xnmctest(:,i+1),Xnmctest(:,i+2),'bs')
    end
    if isempty(Xtruth)==0
        plot(Xntruth(:,i+1),Xntruth(:,i+2),'k*')
    end
    title(['true points and MC: time step = ',num2str(Tstepk1)])
    hold off
end

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
%% GMM HULL
LBtest=-1.5*ones(1,dim);
UBtest=1.5*ones(1,dim);
Xineq=[];
% [Xineq1,~] = GLgn_pts(-1.5*ones(1,dim),1.5*ones(1,dim),5);
Xineq = mvurnd(LBtest,UBtest,20000);

disp(['Running GMM hull'])
GMMHull = GMMFitDataSet(dsX.X,dsX.p);
GMMHull.SetGMM_Hull(15);

indbnd = GMMHull.IsInsideHull(Xineq,1.3);
Xineq = Xineq(~indbnd,:);
figure(26)
GMMHull.plotGMMpointsHUll([1,2],Xineq,1.2,'ro')
hold on
plot3(Xineq(:,1),Xineq(:,2),-ones(size(Xineq,1),1),'b+')
hold off
title('GMM HUll')

%%
XineqTree = mvurnd(LBtest,UBtest,200000);
indbnd = GMMHull.IsInsideHull(XineqTree,1.3);
XineqTree = XineqTree(~indbnd,:);

RigTree=RegInterpFitters('DecisionTreeAdaptiveOutRegion');
RigTree.fit(dsX.X,dsX.p,XineqTree,[],[],GMMHull)
% Xineq = getXineq6D(RigTree,GMMHull,20000,LBtest,UBtest);
% indbnd = GMMHull.IsInsideHull(Xineq,1.3);
% Xineq = Xineq(~indbnd,:);

%% fitting poly to log of probas
% keyboard
% Nm = 3;

Pf=Basis_polyND(dim,4);
lamdim = length(Pf);


%% reguralization points

% keyboard
%% 
% close all
% 
% LB1 = -0.1*ones(dim,1);
% UB1 = 0.1*ones(dim,1);
% 
% LB2 = -1*ones(dim,1);
% UB2 = 1*ones(dim,1);
% 
% % LB = -0.9*ones(dim,1);
% % UB = 0.9*ones(dim,1);
% 
% [Xnew,w]=GH_pts(zeros(dim,1),eye(dim),5);
% Xnew1 = boxShift_working(Xnew,min(Xnew,[],1),max(Xnew,[],1),LB1,UB1);
% Xnew2 = boxShift_working(Xnew,min(Xnew,[],1),max(Xnew,[],1),LB2,UB2);
% 
% % [Xnew1,w]=GLgn_pts(LB1,UB1,5);
% % [Xnew2,w]=GLgn_pts(LB2,UB2,5);
% 
% xpoint= dsX.X(1,:);
% 
% Xnew=[Xnew1;Xnew2];
% Xnew(:,1)=Xnew(:,1);
% dd=xpoint-mean(Xnew,1);
% Xnew = Xnew + repmat(dd,size(Xnew,1),1);
% 

% [Xntruek1,~] = dsX.ApplyAffineTransform_Final2Original(Xnew, zeros(size(Xnew,1),1) );
% 
% % Xntruek=zeros(size(Xntruek1));
% pntruek1=zeros(size(Xnew,1),1);

% for i=1:size(Xnew,1)
%     Xntruek=model.fback(dtkk1,Tk1,Xntruek1(i,:));
%     XnNormk = priorpdfnormk.transForms.trueX2normX(Xntruek(:));
%     pnnormk = priorpdfnormk.func(XnNormk);
%     if isinf(pnnormk)
%         keyboard
%     end
%     pntruek1(i) = priorpdfnormk.transForms.normprob2trueprob(pnnormk);    
% end
% 
% pnew = dsX.ApplyAffineTransformProb_Original2Final(pntruek1);


% [~,ind]=sort(pnew(:),1,'descend');
% pnew=pnew(ind);
% Xnew=Xnew(ind,:);
% 
% indd = pnew>1e-300;
% pnew = pnew(indd);
% Xnew=Xnew(indd,:);

Xnew = dsX.X;
pnew = dsX.p;
% indd = pnew<1e-70;
% pnew(indd) = 1e-70;

for i=0:4
    figure(33+i)
    dsX.PlotPointsProbs3D([i+1,i+2],'ro');
    hold on
    plot3(Xnew(:,i+1),Xnew(:,i+2),pnew,'g+')
    if isempty(Xmctest)==0
        plot(Xnmctest(:,i+1),Xnmctest(:,i+2),'bs')
    end
    if isempty(Xtruth)==0
        plot(Xntruth(:,i+1),Xntruth(:,i+2),'k*')
    end
    title(['true points and MC: time step = ',num2str(Tstepk1)])
    hold off
end

% figure(33)
% hold on
% plot3(Xnew(:,1),Xnew(:,2),pnew,'g+')
% hold off


factconst = max(pnew)/10;
pnfit = pnew/factconst;
beq = log(pnfit);

Aeq=zeros(size(Xnew,1),lamdim);
for ib=1:length(Pf)
    Aeq(:,ib)=evaluate_PolyforLargetSetX(Pf{ib},Xnew);
end

Aineq=zeros(size(Xineq,1),lamdim);
for ib=1:length(Pf)
    Aineq(:,ib)=evaluate_PolyforLargetSetX(Pf{ib},Xineq);
end
bineq = min(beq)*ones(size(Xineq,1),1);
Nineq = size(Xineq,1);

NN = size(Xnew,1);

NNeq=10;
AeqTop=Aeq(1:NNeq,:);
beqTop=beq(1:NNeq);

C=0.03;
% wts=pnew/sum(pnew);
wts=ones(size(pnew));
wts(1:100)=20;

cvx_begin
    variables teq(NN) lam(lamdim) t%  teqTop(Ntop)
    minimize( C*norm(lam,1)+1*norm(teq.*(wts),2) +t ) % +500*norm(teqTop,2)
    subject to
    Aeq*lam==beq+teq;
    Aineq*lam<=bineq+t*ones(Nineq,1);
%     t>=0;
%     AeqTop*lam==beqTop; %+teqTop
cvx_end

% lam = Aeq\beq


lam(1) = lam(1)+log(factconst);
lamsol = lam;

mxentpoly=zeros(1,dim+1);
for i=1:lamdim
    mxentpoly=add_sub_polyND(mxentpoly, scalar_multiply_polyND(lamsol(i),Pf{i}),'add');
end
mxentpoly=simplify_polyND(mxentpoly);

%

pdfnorm.dim =dim;
pdfnorm.mxentpoly_norm = mxentpoly;

LB = -1.2*ones(dim,1);
UB = 1.2*ones(dim,1);
pdfnorm.func=@(x)exppdf_withbounds(x,mxentpoly,LB,UB);
pdfnorm.transForms = dsX.GetTrasnformers();

pdfnorm.info = 'true-0I-hypercube-11';
pdfnorm.pdftype = 'ExpPdf';

pdfnorm.GMMHull = GMMHull;
pdfnorm.LB = LB;
pdfnorm.UB = UB;
pdfnorm.RigTree = RigTree;

% normconst = integratorFuncTrueX_usingpdfnorm(pdfnorm,@(x)constantfunc(x,1),'RegTreeBoxIntegrator');
% 
% pdfnorm.func=@(x)(1/normconst)*exp(evaluate_polyND(mxentpoly_norm,x));

%
Xtest = [mvurnd(-1.2*ones(dim,1),1.2*ones(dim,1),200000);mvurnd(-0.6*ones(dim,1),0.6*ones(dim,1),20000)];

pestXnew = pdfnorm.func(Xnew);
pestXtest = pdfnorm.func(Xtest);
for i=0:4
    figure(33+i)
    dsX.PlotPointsProbs3D([i+1,i+2],'ro');
    hold on
    plot3(Xnew(:,i+1),Xnew(:,i+2),pnew,'g+')
    plot3(Xnew(:,i+1),Xnew(:,i+2),pestXnew,'bs')
    plot3(Xtest(:,i+1),Xtest(:,i+2),pestXtest,'k.')
    if isempty(Xmctest)==0
        plot(Xnmctest(:,i+1),Xnmctest(:,i+2),'bs')
    end
    if isempty(Xtruth)==0
        plot(Xntruth(:,i+1),Xntruth(:,i+2),'k*')
    end
    title(['true points and MC: time step = ',num2str(Tstepk1)])
    hold off
end

% figure(33)
% hold on
% plot3(Xnew(:,1),Xnew(:,2),pdfnorm.func(Xnew),'bs')
% plot3(Xtest(:,1),Xtest(:,2),pdfnorm.func(Xtest),'k.')
% hold off
%%

% keyboard
