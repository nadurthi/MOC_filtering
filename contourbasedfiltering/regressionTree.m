clear 
close all
clc
load('WorstCaseDuffTest.mat')
% X=X(:,[1,2]);
% mX=mX(1:2);
% PX=PX(1:2,1:2);
[mX,PX]=MeanCov(X,probs/sum(probs));
dsX = DataSet(X,probs,'TrueState');
[N,dim] =size(X);


dsX.AddMeanCov_to_OI_Trasform(mX,PX);
dsX.AddHyperCubeTrasform(-1*ones(dim,1),1*ones(dim,1));

indd = dsX.p<1e-70;
dsX.p(indd) = 1e-70;

dsX.SortByProb('descend');

logpn = dsX.getlogProb();
figure
plot3(dsX.X(:,1),dsX.X(:,2),dsX.p,'ro')
% figure
% plot3(dsX.X(:,4),dsX.X(:,5),dsX.p,'ro')
%%
GMMhull= GMMFitDataSet(dsX.X,dsX.p);
% GMM = GMMfitter.FitGMM_kmeans_optimwt(5);
Xineq=[];
% [Xbnd1,~] = GLgn_pts(-1.5*ones(1,dim),1.5*ones(1,dim),4);

SS = mvurnd(-1.5*ones(dim,1),1.5*ones(dim,1),10000);

Xineq = [Xineq;SS];

GMMhull.SetGMM_Hull(15);
indbnd = GMMhull.IsInsideHull(Xineq,1.5);
Xineq = Xineq(~indbnd,:);


close all
GMMhull.plotGMMpointsHUll([1,2],Xineq,1.2,'ro')
hold on
RIF=RegInterpFitters('DecisionTreeAdaptiveOutRegion');
RIF.fit(dsX.X,dsX.p,Xineq,[],[],GMMhull)
RIF.plot(dsX.X,dsX.p,[],[1,2],-1.5*ones(1,dim),1.5*ones(1,dim))
hold on
plot2Dboxes(RIF.method_params.tree)
pest= RIF.evalfit(dsX.X);
plot3(dsX.X(:,1),dsX.X(:,2),dsX.p,'ro')
plot3(dsX.X(:,1),dsX.X(:,2),pest,'b+')

SS = mvurnd(-1.5*ones(dim,1),1.5*ones(dim,1),10000);
pSSest= RIF.evalfit(SS);
plot3(SS(:,1),SS(:,2),pSSest,'gs')

%%

close all
RIF=RegInterpFitters('KnnMean');
RIF.fit(dsX.X,dsX.p,Xineq,[],[])
RIF.plot(dsX.X,dsX.p,Xineq,[1,2],-1.5*ones(1,dim),1.5*ones(1,dim))
%%

close all
RIF=RegInterpFitters('KnnLeastDeg');
RIF.fit(dsX.X,dsX.p,Xineq,[],[])
RIF.plot(dsX.X,dsX.p,Xineq,[1,2],-1.5*ones(1,dim),1.5*ones(1,dim))
%%
close all
Pf=Basis_polyND(dim,1);

RIF=RegInterpFitters('KnnPolyFit');
RIF.fit(dsX.X,dsX.p,Xineq,[],Pf)
RIF.plot(dsX.X,dsX.p,Xineq,[1,2],-1.5*ones(1,dim),1.5*ones(1,dim))
%%

close all
RIF=RegInterpFitters('CARTgmm');
RIF.fit(dsX.X,dsX.p,Xineq,[],[])
RIF.plot(dsX.X,dsX.p,Xineq,[1,2],-1.5*ones(1,dim),1.5*ones(1,dim))

%%
% 
clc
% gprMdl = fitrgp(dsX.X(1:10000,:),log(dsX.p(1:10000,:)),'Basis','linear','FitMethod','exact','PredictMethod','exact') ;
% ypred = resubPredict(gprMdl);

disp('done')
%%

% mp = mpapi(dsX.X(:,[1,2])', log(dsX.p)');
% if(mpval(mp,points)

nc=20;
Z = linkage(dsX.X,'centroid');
T = cluster(Z,'maxclust',nc);
disp('done')

figure
hold on
for n=1:nc
    plot(dsX.X(T==n,4),dsX.X(T==n,5),'bo')
    axis([-2,2,-2,2])
    pause(1)
end
%%
nc=50;
[idx,C] = kmeans(dsX.X,nc);
disp('done')

figure
hold on
for n=1:nc
    plot(dsX.X(idx==n,4),dsX.X(idx==n,5),'bo')
    axis([-2,2,-2,2])
    pause(1)
end

%%
clc
% states={1,2,3,4,5,6};
states={'x','y','z','vx','vy','vz'};
states={'1','2','3','4','5','6'};
tree = fitrtree(dsX.X,dsX.p,'MinParentSize',10,'MaxNumSplits',10000,'MinLeafSize',10,...
                'PredictorNames',states,'ResponseName','probs');
            
6+6
%%
tree.Children
tree.CutPoint
tree.CutPredictor
tree.IsBranchNode
tree.NodeMean
tree.Parent
%%
close all
pp = predict(tree,dsX.X);
figure
plot3(dsX.X(:,1),dsX.X(:,2),pp,'bo')
hold on
plot3(dsX.X(:,1),dsX.X(:,2),(dsX.p),'ro')

figure
plot3(dsX.X(:,4),dsX.X(:,5),pp,'bo')
hold on
plot3(dsX.X(:,4),dsX.X(:,5),(dsX.p),'ro')
%% plot the cuts
% min,max,#node,type{node,leaf} ...
boxes={-1.2*ones(1,dim),1.2*ones(1,dim),1};
tree=RIF.method_params.tree;
for nc = 1:length(tree.CutPredictor)
    st = tree.CutPredictor{nc};
    if strcmp(tree.CutPredictor{nc},'')==false % node has a branch
        cp = tree.CutPoint(nc);
        stind = getindex2cell(states,st);
        for r=1:nc+10
            if boxes{r,3}==nc
                break
            end
        end
        
        boxmin = boxes{r,1};
        boxmax = boxes{r,2};
        
        newminL=boxmin;
        newmaxL = boxmax;
        newmaxL(stind) = cp;
        
        newminR=boxmin;
        newmaxR = boxmax;
        newminR(stind) = cp;
        
        l = size(boxes,1);
        boxes{l+1,1}=newminL;
        boxes{l+1,2}=newmaxL;
        boxes{l+1,3}=tree.Children(nc,1);
        
        boxes{l+2,1}=newminR;
        boxes{l+2,2}=newmaxR;
        boxes{l+2,3}=tree.Children(nc,2);
        
        boxes(r,:)=[];
        
    else  % node has no brachs or cildern ... hence it is a leaf
        
    end
    
end
disp('done')
%%
close all
tree = RIF.method_params.tree;
boxes=getTree2Boxes(tree);
D=[]
for r=1:size(boxes,1)
    mnb = boxes{r,1};
    mxb = boxes{r,2};
    dr=getALLcorners(mnb,mxb);
    D=[D;max(abs(mxb-mnb))];
    for s=1:1
        figure(s)
%         plot3(dsX.X(:,s),dsX.X(:,s+1),pp(:),'bo')
%         hold on
        plot3(dr(:,s),dr(:,s+1),-1*ones(size(dr,1),1),'r','linewidth',2)
        axis([-2,2,-2,2])
        title(num2str(r))
        hold off
    end
    pause(0.2)
    
end








