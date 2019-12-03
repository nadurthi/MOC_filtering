load('6Dtestcase_backprop.mat')
dim=6;
dtkk1 = time.dt;

% inmain: model, priorfullnormpdf, X,   probs,   4,  k,       time.Tvec(k), time.dt, Xmctest, Xtruth(k,:)
% infunc: model, priorfullnormpdf, Xk1, probsk1, Nm, Tstepk1, Tk1,          dtkk1,   Xmctest, Xtruth
%%
    Xmctestnorm = priorfullnormpdf.transForms.trueX2normX(Xmctest);
    pMCnotm = priorfullnormpdf.transForms.trueprob2normprob(pMC) ;
    
%     pMCestnorm = exp(evaluate_polyND(PnD,Xmctestnorm)); 
    pMCestnorm = priorfullnormpdf.func(Xmctestnorm);
    
    close all
    figure
    plot3(Xmctestnorm(:,1),Xmctestnorm(:,2),pMCnotm,'r+',Xmctestnorm(:,1),Xmctestnorm(:,2),pMCestnorm,'bo')
    figure
    plot3(Xmctestnorm(:,2),Xmctestnorm(:,3),pMCnotm,'r+',Xmctestnorm(:,2),Xmctestnorm(:,3),pMCestnorm,'bo')
    figure
    plot3(Xmctestnorm(:,3),Xmctestnorm(:,4),pMCnotm,'r+',Xmctestnorm(:,3),Xmctestnorm(:,4),pMCestnorm,'bo')
    figure
    plot3(Xmctestnorm(:,4),Xmctestnorm(:,5),pMCnotm,'r+',Xmctestnorm(:,4),Xmctestnorm(:,5),pMCestnorm,'bo')
    figure
    plot3(Xmctestnorm(:,5),Xmctestnorm(:,6),pMCnotm,'r+',Xmctestnorm(:,5),Xmctestnorm(:,6),pMCestnorm,'bo')
    
    
%% test if it blows up
hh=0.9;
LB=-hh*ones(dim,1);
UB=hh*ones(dim,1);

PnD = priorfullnormpdf.mxentpoly_norm;
% first get the test points
XcheckA = mvurnd(LB,UB,dim*1000);
pcheckA = priorfullnormpdf.func(XcheckA);

close all
figure
hist(pcheckA,1000)
figure
hist(pMCnotm,1000)
%%

[Xqi,pts1Dqi,interpPoly1Dqi]=getsparsePts(dim+5,dim);

pcheckB = priorfullnormpdf.func(Xqi);
close all
figure
hist(pcheckB,1000)

%%
clc
clear all
close all

load('6Dtestcase_backprop.mat')
dim=6;
dtkk1 = time.dt;


kk=3;
Tstepk1=kk;
Tk1=time.Tvec(kk);

[X,probs]=propagate_character(X,probs,time.dt,time.Tvec(kk),model);
[aa,bb]=MeanCov(X,probs/sum(probs));

Xmctest = zeros(size(XMC,1),model.fn);
for ii=1:size(XMC,1)
    Xmctest(ii,:) = XMC(ii,:,kk);
end

[N,dim] =size(Xmctest);

[mXk1,PXk1]=MeanCov(Xmctest,1/N*ones(N,1));


% Xmctestnormk=Xmctestnorm;
% pMCnotmk=pMCnotm;
%%
dsX = DataSet(X,probs,'TrueState');

PP = diag(diag(PXk1));
dsX.AddMeanCov_to_OI_Trasform(mXk1,1^2*PXk1);
% dsX.AddHyperCubeTrasform(-1*ones(dim,1),1*ones(dim,1));

% plot3(dsX.X(:,1),dsX.X(:,2),dsX.p,'ro')

% keyboard

% indd = dsX.p<1e-100;
% dsX.p(indd) = 1e-100;

dsX.SortByProb('descend');

logpn = dsX.getlogProb();

if isempty(Xmctest)==0
    [Xnmctest,pMCn] = dsX.ApplyAffineTransform_Original2Final(Xmctest, pMC );
end
Xtruth=[];
if isempty(Xtruth)==0
    [Xntruth,~] =  dsX.ApplyAffineTransform_Original2Final(Xtruth, zeros(size(Xtruth,1),1) );
end


%% plottinmg
close all
for i=0:4
    figure(80+i)
    dsX.PlotPointsProbs3D([i+1,i+2],'ro');
    hold on
    if isempty(Xmctest)==0
        plot3(Xnmctest(:,i+1),Xnmctest(:,i+2),pMCn,'bs')
    end
    if isempty(Xtruth)==0
        plot(Xntruth(:,i+1),Xntruth(:,i+2),'k*')
    end
    title(['true points and MC: time step = ',num2str(Tstepk1)])
    hold off
end


%% Form the hyepercube and go get back the fucntion valuess
% startpdfdata = pdftraj{Tstepk1,1};
transForms_atk = priorfullnormpdf.transForms;
normpdf_atk = priorfullnormpdf.func;
normpdfexp_atk = priorfullnormpdf.expfunc;
% normpdfexp_atk = @(x)evaluate_polyND(priorfullnormpdf.mxentpoly_norm,x);
transForms_atk1 = dsX.GetTrasnformers();

LB=-1*ones(dim,1);
UB=1*ones(dim,1);

% first get the test points
Xcheckk1 = mvurnd(1*LB,1*UB,dim*500);
probscheck1 = zeros(size(Xcheckk1,1),1);
logprobscheck1 = zeros(size(Xcheckk1,1),1);

xprev=zeros(size(Xcheckk1));
    xnormprev=zeros(size(Xcheckk1));
    
for i=1:size(Xcheckk1,1)
    xnorm = Xcheckk1(i,:);
    [probscheck1(i),logprobscheck1(i),xprev(i,:),xnormprev(i,:)] = getbacknormprobs(model,Tk1,dtkk1,normpdf_atk,normpdfexp_atk,xnorm,transForms_atk1,transForms_atk);
end
% logprobscheck1=log(probscheck1+1e-12);
close all
figure(71)
plot3(Xcheckk1(:,1),Xcheckk1(:,2),probscheck1,'bs')
hold on
plot3(Xnmctest(:,1),Xnmctest(:,2),pMCn,'ro')
hold off


for i=0:4
figure
plot(xnormprev(:,i+1),xnormprev(:,i+2),'bs',Xmctestnormk(:,i+1),Xmctestnormk(:,i+2),'ro')
end

 LB=-1.1*ones(dim,1);
    UB=1.1*ones(dim,1);
    insidept = isinbox(xnormprev,LB,UB);
    sum(insidept)
%%
Xall=[];
Fall=[];
logFall=[];
preverr=1e10;
prevans=[];
for qi=4
    [Xqi,pts1Dqi,interpPoly1Dqi]=getsparsePts_scaled(dim+qi,dim,1);
%     Xqi = 0.6*Xqi;

    if isempty(Xall)
        %         Xall = Xqi;
        Xrest = Xqi;
    else
        Xrest = setdiff(Xqi,Xall,'rows');
    end
    %     Xk1qi = zeros(size(Xqi));
    probsk1qi = zeros(size(Xrest,1),1);
    logprobsk1qi = zeros(size(Xrest,1),1);
    xprev=zeros(size(Xrest));
    xnormprev=zeros(size(Xrest));
    
    for i=1:size(Xrest,1)
        [probsk1qi(i),logprobsk1qi(i),xprev(i,:),xnormprev(i,:)] = getbacknormprobs(model,Tk1,dtkk1,normpdf_atk,normpdfexp_atk,Xrest(i,:)',transForms_atk1,transForms_atk);
    end
    Xall = vertcat(Xall,Xrest);
    Fall = vertcat(Fall,probsk1qi);
    logFall = vertcat(logFall,logprobsk1qi);
    %     keyboard
    
    LB=-1.1*ones(dim,1);
    UB=1.1*ones(dim,1);
    insidept = isinbox(xnormprev,LB,UB);
    
    FuncTablelog = [Xall,logFall];
    
    % figure
    % plot3(X(:,1),X(:,2),FuncTable(:,3),'bo')
    
    PnD=sparseProductInterpPoly(dim+qi,dim,pts1Dqi,interpPoly1Dqi,FuncTablelog,'NestedClenshawCurtis');
    
    Festlog=evaluate_polyND(PnD,Xall);
    [Festlog,FuncTablelog(:,end)]
    
    Flogcheck = evaluate_polyND(PnD,Xcheckk1);
    Fcheck = exp(Flogcheck);
    max(abs((logprobscheck1-Flogcheck)))
    ind = abs(logprobscheck1)>1e-6;
    err=max( abs((logprobscheck1(ind)-Flogcheck(ind))./logprobscheck1(ind)) );
    %     keyboard
    
    if 100*err < 1
        %         prevans = PnD;
        break
    end
    %     if err > preverr
    %         prevans = PnD;
    %         break
    %     end
    %     preverr = err;
    %     prevans = PnD;
    
end
% PnD=prevans;
[sum(insidept),length(insidept)]

close all
for i=0:4
    figure
    plot(xnormprev(:,i+1),xnormprev(:,i+2),'bs',Xmctestnormk(:,i+1),Xmctestnormk(:,i+2),'ro')
end
% Xmctestnormk=Xmctestnorm;
% pMCnotmk=pMCnotm;

%% truth test
logprobest_truth=evaluate_polyND(PnD,Xnmctest);
probest_truth = exp(logprobest_truth);
[probest_truth,pMCn]
figure(89)
plot3(Xnmctest(:,1),Xnmctest(:,2),probest_truth,'bo',Xnmctest(:,1),Xnmctest(:,2),pMCn,'r+')

% keyboard
%% debug plots
if false
    [xx,yy]=meshgrid(linspace(-1,1,25));
    Fplot12 = zeros(size(xx));
    Fplot23 = zeros(size(xx));
    Fplot34 = zeros(size(xx));
    Fplot45 = zeros(size(xx));
    Fplot56 = zeros(size(xx));
    
    a12 = cell(size(xx,1),1);
    a23 = cell(size(xx,1),1);
    a34 = cell(size(xx,1),1);
    a45 = cell(size(xx,1),1);
    a56 = cell(size(xx,1),1);
    
    parfor i=1:size(xx,1)
        a12{i}=zeros(size(xx,2),1);
        a23{i}=zeros(size(xx,2),1);
        a34{i}=zeros(size(xx,2),1);
        a45{i}=zeros(size(xx,2),1);
        a56{i}=zeros(size(xx,2),1);
        
        for j=1:size(xx,2)
            [i,j]
            %         a12{i}(j) = marginalize_exp_pdf_method2([xx(i,j),yy(i,j)],[1,2],@(g)exp(evaluate_polyND(PnD,g)),dim);
            
            a12{i}(j) = marginalize_exp_pdf_method2([xx(i,j),yy(i,j)],[1,2],@(g)exp(evaluate_polyND(PnD,g)),dim);
            %         a23{i}(j) = marginalize_exp_pdf_method2([xx(i,j),yy(i,j)],[2,3],@(g)exp(evaluate_polyND(PnD,g)),dim);
            %         a34{i}(j) = marginalize_exp_pdf_method2([xx(i,j),yy(i,j)],[3,4],@(g)exp(evaluate_polyND(PnD,g)),dim);
            %         a45{i}(j) = marginalize_exp_pdf_method2([xx(i,j),yy(i,j)],[4,5],@(g)exp(evaluate_polyND(PnD,g)),dim);
            %         a56{i}(j) = marginalize_exp_pdf_method2([xx(i,j),yy(i,j)],[5,6],@(g)exp(evaluate_polyND(PnD,g)),dim);
            
            %         Fplot12(i,j) = evaluate_polyND(PnD,[xx(i,j),yy(i,j),0,0,0,0]);
            %         Fplot23(i,j) = evaluate_polyND(PnD,[0,xx(i,j),yy(i,j),0,0,0]);
            %         Fplot34(i,j) = evaluate_polyND(PnD,[0,0,xx(i,j),yy(i,j),0,0]);
            %         Fplot45(i,j) = evaluate_polyND(PnD,[0,0,0,xx(i,j),yy(i,j),0]);
            %         Fplot56(i,j) = evaluate_polyND(PnD,[0,0,0,0,xx(i,j),yy(i,j)]);
        end
    end
    for i=1:size(xx,1)
        Fplot12(i,:)= a12{i};
        Fplot23(i,:)= a23{i};
        Fplot34(i,:)= a34{i};
        Fplot45(i,:)= a45{i};
        Fplot56(i,:)= a56{i};
    end
    Fplot12=exp(Fplot12)-0;
    Fplot23=exp(Fplot23)-0;
    Fplot34=exp(Fplot34)-0;
    Fplot45=exp(Fplot45)-0;
    Fplot56=exp(Fplot56)-0;
    
    figure(90)
    surf(xx,yy,Fplot12,'FaceColor','r','EdgeColor','none','FaceAlpha',0.5);
    camlight right; lighting phong
    alpha 0.5
    figure(91)
    contour(xx,yy,Fplot12,15);
    
    
    figure(92)
    surf(xx,yy,Fplot23,'FaceColor','r','EdgeColor','none','FaceAlpha',0.5);
    camlight right; lighting phong
    alpha 0.5
    figure(93)
    contour(xx,yy,Fplot23,15);
    
    figure(94)
    surf(xx,yy,Fplot34,'FaceColor','r','EdgeColor','none','FaceAlpha',0.5);
    camlight right; lighting phong
    alpha 0.5
    figure(95)
    contour(xx,yy,Fplot34,15);
    
    figure(96)
    surf(xx,yy,Fplot45,'FaceColor','r','EdgeColor','none','FaceAlpha',0.5);
    camlight right; lighting phong
    alpha 0.5
    figure(97)
    contour(xx,yy,Fplot45,15);
    
    figure(98)
    surf(xx,yy,Fplot56,'FaceColor','r','EdgeColor','none','FaceAlpha',0.5);
    camlight right; lighting phong
    alpha 0.5
    figure(99)
    contour(xx,yy,Fplot56,15);
    
end
%% normalize
xl = -1*ones(1,dim);
xu =  1*ones(1,dim);
vv = prod(xu-xl);
for Ninteg=7
    [xint,wint] = GLeg_pts(Ninteg*ones(1,dim), xl, xu);
    pp=exp(evaluate_polyND(PnD,xint));
    cc=vv*wint(:)'*pp(:);
end
disp(['norm constant is = ',num2str(cc)])
keyboard
% c0 = get_coeff_NDpoly(PnD,zeros(1,dim));
% PnD = update_or_insert_coeff_NDpoly(PnD,zeros(1,dim),c0-log(cc));
%%
pdfnorm.dim =dim;
pdfnorm.mxentpoly_norm = PnD;

LB = -1.1*ones(dim,1);
UB = 1.1*ones(dim,1);
pdfnorm.func=@(x)exppdf_withbounds(x,PnD,LB,UB);
pdfnorm.expfunc=@(x)expofaexppdf_withbounds(x,PnD,LB,UB);

pdfnorm.transForms = dsX.GetTrasnformers();

pdfnorm.info = 'true-50I';
pdfnorm.pdftype = 'ExpPdf';

% pdfnorm.GMMHull = GMMHull;
pdfnorm.LB = LB;
pdfnorm.UB = UB;

% [Xqi,pts1Dqi,interpPoly1Dqi]=getsparsePts(dim+3,dim);
[Xqi,~] = GH_points(zeros(dim,1),(1/3)^2*eye(dim),4);
% mm=max(max(Xqi));
% Xqi = 1.1*Xqi/mm;
normprobs = pdfnorm.func(Xqi);
Xnew=pdfnorm.transForms.normX2trueX(Xqi);
probsnew=pdfnorm.transForms.normprob2trueprob(normprobs);


figure(100)
plot3(Xnew(:,1),Xnew(:,2),probsnew,'bo',Xk1(:,1),Xk1(:,2),probsk1,'r+')


% pdfnorm.RigTree = RigTree;

% normconst = integratorFuncTrueX_usingpdfnorm(pdfnorm,@(x)constantfunc(x,1),'RegTreeBoxIntegrator');
%
% pdfnorm.func=@(x)(1/normconst)*exp(evaluate_polyND(mxentpoly_norm,x));

%
