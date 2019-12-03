function [pdfnorm,Xnew,probsnew] = get_interp_pdf_0I_6D_gobackget_coe(model,priorfullnormpdf,Xk1,probsk1,Nm,Tstepk1,Tk1,dtkk1,Xmctest,Xtruth)
% known prior pdf priorpdfnormk at time k
% compute pdfnorm or pdfnormk1 at time k+1
% Tk1 is time at k+1,
% dtkk1 is deltat from k to k+1

%%
[N,dim] =size(Xk1);

[mXk1,PXk1]=MeanCov(Xk1,probsk1/sum(probsk1));


dsX = DataSet(Xk1,probsk1,'TrueState');


dsX.AddMeanCov_to_OI_Trasform(mXk1,3^2*PXk1);
% dsX.AddHyperCubeTrasform(-1*ones(dim,1),1*ones(dim,1));

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
    figure(80+i)
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


%% Form the hyepercube and go get back the fucntion valuess
% startpdfdata = pdftraj{Tstepk1,1};
transForms_atk = priorfullnormpdf.transForms;
normpdf_atk = priorfullnormpdf.func;
normpdfexp_atk = priorfullnormpdf.expfunc;
transForms_atk1 = dsX.GetTrasnformers();

LB=-ones(dim,1);
UB=ones(dim,1);

% first get the test points
Xcheckk1 = mvurnd(0.8*LB,0.8*UB,dim*5000);
probscheck1 = zeros(size(Xcheckk1,1),1);
logprobscheck1 = zeros(size(Xcheckk1,1),1);

for i=1:size(Xcheckk1,1)
    xnorm = Xcheckk1(i,:);
    %     [xk1,~]=transForms_atk1.normX2trueX(x) ;
    %
    %     xk = model.fback(dtkk1,Tk1,xk1);
    %     xknorm = transForms_atk.trueX2normX(xk);
    %     pknorm = normpdf_atk(xknorm);
    %     pk = transForms_atk.normprob2trueprob(pknorm);
    %     pk1 = pk; % as sat dynamics do not diverge
    %     probscheck1(i)=transForms_atk1.trueprob2normprob(pk1);
    [probscheck1(i),logprobscheck1(i)] = getbacknormprobs(model,Tk1,dtkk1,normpdf_atk,normpdfexp_atk,xnorm,transForms_atk1,transForms_atk);
end
% logprobscheck1=log(probscheck1+1e-12);
figure(71)
plot3(Xcheckk1(:,1),Xcheckk1(:,2),probscheck1,'bo')

Xall=[];
Fall=[];
logFall=[];
preverr=1e10;
prevans=[];
for qi=4
    [Xqi,pts1Dqi,interpPoly1Dqi]=getsparsePts(dim+qi,dim);
    
    if isempty(Xall)
        %         Xall = Xqi;
        Xrest = Xqi;
    else
        Xrest = setdiff(Xqi,Xall,'rows');
    end
    %     Xk1qi = zeros(size(Xqi));
    probsk1qi = zeros(size(Xrest,1),1);
    logprobsk1qi = zeros(size(Xrest,1),1);
    xkvec = zeros(size(Xrest,1),dim);
    xknormvec = zeros(size(Xrest,1),dim);
    
    for i=1:size(Xrest,1)
        [probsk1qi(i),logprobsk1qi(i),xkvec(i,:),xknormvec(i,:)] = getbacknormprobs(model,Tk1,dtkk1,normpdf_atk,normpdfexp_atk,Xrest(i,:)',transForms_atk1,transForms_atk);
    end
    Xall = vertcat(Xall,Xrest);
    Fall = vertcat(Fall,probsk1qi);
    logFall = vertcat(logFall,logprobsk1qi);
    insidept = isinbox(xknormvec,LB,UB);
    
    FuncTablelog = [Xall,logFall];
%     ind=FuncTablelog(:,3)<-300;
%     FuncTablelog(ind,3)=-30;
    % figure
    % plot3(X(:,1),X(:,2),FuncTable(:,3),'bo')
    
    PnD=sparseProductInterpPoly(dim+qi,dim,pts1Dqi,interpPoly1Dqi,FuncTablelog,'NestedClenshawCurtis');
    
    Festlog=evaluate_polyND(PnD,Xall);
    [Festlog,FuncTablelog(:,end)]
    
    Flogcheck = evaluate_polyND(PnD,Xcheckk1);
    Fcheck = exp(Flogcheck);
    max(abs((logprobscheck1-Flogcheck)))
    max(abs((exp(logprobscheck1)-exp(Flogcheck))))
    ind = abs(logprobscheck1)>1e-6;
    err=max( abs((logprobscheck1(ind)-Flogcheck(ind))./logprobscheck1(ind)) );
    err2=max( abs((exp(logprobscheck1(ind))-exp(Flogcheck(ind)))./exp(logprobscheck1(ind))) );
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
Fall(logprobsk1qi>5)
Xall(logprobsk1qi>5,:)
xknormvec(logprobsk1qi>5,:)
[sum(insidept),length(insidept)]
%% truth test
logprobest_truth=evaluate_polyND(PnD,dsX.X);
% f=evaluate_polyND_3(PnD,dsX.X)
probest_truth = exp(logprobest_truth);
[probest_truth,dsX.p]
figure(89)
plot3(dsX.X(:,1),dsX.X(:,2),probest_truth,'bo',dsX.X(:,1),dsX.X(:,2),dsX.p,'r+')

Xcheckk2 = mvurnd(-1.0*ones(1,dim),1.0*ones(1,dim),dim*10000);
probest_check2 = exp(evaluate_polyND(PnD,Xcheckk2));
figure(69)
plot3(Xcheckk2(:,1),Xcheckk2(:,2),probest_check2,'ro')

keyboard
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
% xl = -1*ones(1,dim);
% xu =  1*ones(1,dim);
% vv = prod(xu-xl);
% for Ninteg=7
%     [xint,wint] = GLeg_pts(Ninteg*ones(1,dim), xl, xu);
%     pp=exp(evaluate_polyND(PnD,xint));
%     cc=vv*wint(:)'*pp(:);
% end
% disp(['norm constant is = ',num2str(cc)])


% keyboard


% c0 = get_coeff_NDpoly(PnD,zeros(1,dim));
% PnD = update_or_insert_coeff_NDpoly(PnD,zeros(1,dim),c0-log(cc));
%%
pdfnorm.dim =dim;
pdfnorm.mxentpoly_norm = PnD;

LB = -1.0*ones(dim,1);
UB = 1.0*ones(dim,1);
pdfnorm.func=@(x)exppdf_withbounds(x,PnD,LB,UB);
pdfnorm.expfunc=@(x)expofaexppdf_withbounds(x,PnD,LB,UB);

pdfnorm.transForms = dsX.GetTrasnformers();

pdfnorm.info = 'true-50I';
pdfnorm.pdftype = 'ExpPdf';

% pdfnorm.GMMHull = GMMHull;
pdfnorm.LB = LB;
pdfnorm.UB = UB;

% [Xqi,pts1Dqi,interpPoly1Dqi]=getsparsePts(dim+3,dim);
[Xqi,~] = GH_points(zeros(dim,1),(0.3)^2*eye(dim),4);
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
