clc
close all
clear
% polar to cartesian conversion


mur=30;
sigr=2;

muth=pi/2;
sigth=20*pi/180;

xf0 = [mur;muth];
Pf0 = [sigr^2,0;0,sigth^2];

dim = 2;


f=@(x0)[x0(1)*cos(x0(2));x0(1)*sin(x0(2))];
fv=@(x0)[x0(:,1).*cos(x0(:,2)),x0(:,1).*sin(x0(:,2))];
fvback=@(x1)[sqrt(sum(x1.^2,2)),atan(x1(:,2)./x1(:,1))];

detJacfv = @(x0)x0(:,1);

x0mc=mvnrnd(xf0',Pf0,10000);
p0mc=mvnpdf(x0mc,xf0',Pf0);

x1mc = fv(x0mc);
p1mc = p0mc./detJacfv(x0mc);

figure(1)
plot3(x1mc(:,1),x1mc(:,2),p1mc,'r+')
title('x1mc - p1mc')

pdf0.Xrep = mvnrnd(xf0',Pf0,1000);
pdf0.prep = mvnpdf(pdf0.Xrep,xf0',Pf0);


dsX0rep = DataSet(pdf0.Xrep,pdf0.prep,'TrueState');
dsX0rep.AddMeanCov_to_OI_Trasform(xf0(:),Pf0);
% dsX.AddHyperCubeTrasform(-0.9*ones(dim,1),0.9*ones(dim,1));
transForms = dsX0rep.GetTrasnformers();

pdf0.normeval=@(x)mvnpdf(x,zeros(1,dim),eye(dim,dim));
pdf0.transForms = dsX0rep.GetTrasnformers();

pdf0.info = 'true-0I';
pdf0.pdftype = 'ExpPdf';

% pdfnorm.GMMHull = GMMHull;
pdf0.LB = -3*ones(1,dim);
pdf0.UB = 3*ones(1,dim);

figure(2)
subplot(1,2,1)
plot3(pdf0.Xrep(:,1),pdf0.Xrep(:,2),pdf0.prep,'r+')
title('pdf0.Xrep')
subplot(1,2,2)
plot3(dsX0rep.X(:,1),dsX0rep.X(:,2),dsX0rep.p,'r+')
title('dsX0rep.X')

figure(3)
x0mc_norm=pdf0.transForms.trueX2normX(x0mc);
p0mc_norm=pdf0.transForms.trueprob2normprob(p0mc);
subplot(1,2,1)
plot3(x0mc(:,1),x0mc(:,2),p0mc,'r+')
title('x0mc')
subplot(1,2,2)
plot3(x0mc_norm(:,1),x0mc_norm(:,2),p0mc_norm,'r+')
title('x0mc norm - p0mc norm')




%%

[N,dim] =size(pdf0.Xrep);

Xrep1=fv(pdf0.Xrep);
prep1=(pdf0.prep)./detJacfv(pdf0.Xrep);
figure(4)
subplot(1,2,1)
plot3(Xrep1(:,1),Xrep1(:,2),prep1,'r+')
title('Xrep1 (front prop)')
subplot(1,2,2)
plot3(x1mc(:,1),x1mc(:,2),p1mc,'r+')
title('x1mc')



[mX1,PX1]=MeanCov(Xrep1,prep1/sum(prep1));


dsX1 = DataSet(Xrep1,prep1,'TrueState');


% mX1=Xk1(b,:)';
dsX1.AddMeanCov_to_OI_Trasform(mX1,3^2*PX1);
% dsX.AddHyperCubeTrasform(-0.9*ones(dim,1),0.9*ones(dim,1));

% plot3(dsX.X(:,1),dsX.X(:,2),dsX.p,'ro')

% keyboard

% indd = dsX.p<1e-100;
% dsX.p(indd) = 1e-100;


%% Form the hyepercube and go get back the fucntion values for prior pdf


LB=-2*ones(dim,1);
UB=2*ones(dim,1);

X1norm_concat=[];
P1norm_concat=[];
cnt=0;
while 1
    % first get the test points
    if isempty(X1norm_concat)==1
        %             X1norm_newsamp = mvurnd(1*LB,1*UB,dim*100);
        X1norm_newsamp=[dsX1.X;mvurnd_expanded(1*LB,1*UB,100,0)];
        %             P1norm_concat=dsX.p
    else
        % tree based adaptive sampling
        
        nb=size(boxes,1);
        nids=discretesample(errbox, dim*10);
        
        X1norm_newsamp=[];
        for i=1:nb
            ns = sum(nids==i);
            lb = boxes{i,1};
            ub = boxes{i,2};
            X1norm_newsamp=[X1norm_newsamp;mvurnd_expanded(lb,ub,2*ns,0)];
        end
        
        MM=Mbox;
        iig=MM>prctile(MM,65);
        MM(iig)=MM(iig)*2;
        nids=pmf_sampler(MM,dim*10);
        for i=1:nb
            ns = sum(nids==i);
            lb = boxes{i,1};
            ub = boxes{i,2};
            X1norm_newsamp=[X1norm_newsamp;mvurnd_expanded(lb,ub,2*ns,0)];
        end
        
        nids=pmf_sampler(Assratio,dim*10);
        for i=1:nb
            ns = sum(nids==i);
            lb = boxes{i,1};
            ub = boxes{i,2};
            X1norm_newsamp=[X1norm_newsamp;mvurnd_expanded(lb,ub,2*ns,0)];
        end
    end
    Nm_newsamp=size(X1norm_newsamp,1);
    [X1true_newsamp,~] = dsX1.ApplyAffineTransform_Final2Original(X1norm_newsamp,ones(Nm_newsamp,1));
    X0true_newsamp = fvback(X1true_newsamp);
    
    X0norm_newsamp=pdf0.transForms.trueX2normX(X0true_newsamp);
    P0norm_newsamp = pdf0.normeval(X0norm_newsamp);
    
    P0true_newsamp=pdf0.transForms.normprob2trueprob(P0norm_newsamp);
    % for now we take trace of dynamics is 0
    P1true_newsamp = P0true_newsamp./detJacfv(X0true_newsamp);
    
    P1norm_newsamp = dsX1.ApplyAffineTransformProb_Original2Final(P1true_newsamp);
    
    figure(10)
    plot3(X0true_newsamp(:,1),X0true_newsamp(:,2),P0true_newsamp,'bo')
    
    figure(11)
    plot3(X1true_newsamp(:,1),X1true_newsamp(:,2),P1norm_newsamp,'bo')
    
    
    X1norm_concat=[X1norm_concat;X1norm_newsamp];
    P1norm_concat=[P1norm_concat;P1norm_newsamp];
    
    Nconcat = size(X1norm_concat,1);
    
    Ntree=10;
    % [X1norm_concat,P1norm_concat];
    EE=cell(Ntree,2);


    for jj=1:Ntree
        cv = cvpartition(size(X1norm_concat,1),'HoldOut',0.3);
        idx = cv.test;
        
        Xtrain=X1norm_concat(~idx,:);
        ytrain=P1norm_concat(~idx);
        Xtest=X1norm_concat(idx,:);
        ytest=P1norm_concat(idx);
        
        % Now build the tree
        tree = fitrtree(Xtrain,ytrain,'MinParentSize',10,'MaxNumSplits',5000,'MinLeafSize',10);
        ypred=tree.predict(Xtest);
        errtest=max(100*abs(ypred-ytest)./(ytest+1));
        EE{jj,1}=tree;
        EE{jj,2}=errtest;
    end

    EE=sortrows(EE,2);
    tree=EE{1,1};
    errtest=EE{1,2};

    
    boxes=getTree2Boxes(tree,LB,UB);
    
    boxes = sortrows(boxes,4,'descend');
    


    figure(71)
    plot3(X1norm_concat(:,1),X1norm_concat(:,2),P1norm_concat,'bo')
    hold on
    plot3Dtreepatches(boxes)
    hold off

    
    
    errbox=zeros(size(boxes,1),1);
    Mbox=zeros(size(boxes,1),1);
    Assratio=zeros(size(boxes,1),1);

    
    
    for j=1:size(boxes,1)
        lb = boxes{j,1};
        ub = boxes{j,2};
        d=ub-lb;
        Assratio(j) = abs(max(d)/min(d));
        if isnan(Assratio(j)) || isinf(Assratio(j))
            Assratio(j)=1e3;
        end
        if Assratio(j)==0
            Assratio(j)=1;
        end
        if Assratio(j)>1e10
            Assratio(j)=1e9;
        end
        indyy = sum((X1norm_concat> lb(:)') & (X1norm_concat <= ub(:)'),2)==dim;
        a = P1norm_concat(indyy);
        m = boxes{j,4};
        if isempty(a)==1
            errbox(j)=0;
        else
            errbox(j) = mean(abs(a-repmat(m,length(a),1))/(m+1));
        end
        Mbox(j) = m;
    end
    
    

    s=prctile(Mbox,65);
    inds = Mbox>s;
    errbox(inds)=errbox(inds)*3;
    if sum(errbox)~=0
        errbox = errbox/sum(errbox);
    end
    if sum(Assratio)~=0
        Assratio = Assratio/sum(Assratio);
    end
    %         Assratio = Assratio/sum(Assratio);
    
    [~,sortind]=sort(errbox);
    errbox=errbox(sortind);
    Mbox=Mbox(sortind);
    boxes=boxes(sortind,:);
    Assratio=Assratio(sortind,:);
    s=prctile(Assratio,65);
    inds = Assratio>s;
    Assratio(inds)=Assratio(inds)*10;
    if sum(Assratio)~=0
        Assratio = Assratio/sum(Assratio);
    end


    if (errtest<25 && cnt>5) || cnt>15
        break
    end
    cnt=cnt+1;
    [cnt,Nconcat]
    
end

X1norm_concat_prior = X1norm_concat;
P1norm_concat_prior = P1norm_concat;

pdf1.normeval=@(x)treepredict_withbnds(x,tree,LB,UB);
pdf1.tree=tree;
pdf1.boxes=boxes;
pdf1.dim =dim;
pdf1.transForms = dsX1.GetTrasnformers();
pdf1.info = 'prior-pdf-tree';
pdf1.pdftype = 'tree';
pdf1.LB = LB;
pdf1.UB = UB;
% 
% nb=size(boxes,1);
% nids=pmf_sampler(Mbox/sum(Mbox),dim*50);
% XX=[];
% for i=1:nb
%     ns = sum(nids==i);
%     lb = boxes{i,1};
%     ub = boxes{i,2};
%     XX=[XX;mvurnd_expanded(lb,ub,10*ns,0)];
% end
% 
% pp = pdf1.normeval(XX);
% [pdf1.Xrep,pdf1.prep] = dsX.ApplyAffineTransform_Final2Original(XX,pp);
% 
% 
% 
% 
% 
% 
% 
