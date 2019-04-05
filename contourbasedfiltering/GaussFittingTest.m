clear
clc
load('GaussFitting.mat')
%%

[N,dim] =size(Xr);


[mX,PX]=MeanCov(Xr,probs/sum(probs));
 
dsX = DataSet(Xr,probs,'TrueState');


mquad=mX(:);

dsX.AddMeanCov_to_OI_Trasform(mX,PX);
dsX.AddHyperCubeTrasform(-1*ones(dim,1),1*ones(dim,1));

indd = dsX.p<1e-70;
dsX.p(indd) = 1e-70;

dsX.SortByProb('descend');

p = dsX.p;
logpn = dsX.getlogProb();

X=dsX.X;

figure
plot3(X(:,1),X(:,2),dsX.p,'ro')
%%
Xineq=[];
[Xbnd1,~] = GLgn_pts(-1.5*ones(1,dim),1.5*ones(1,dim),8);
Xbnd2=[];
Xineq = [Xbnd1;Xbnd2];


GMMHull = GMMFitDataSet(dsX.X,dsX.p);
GMMHull.SetGMM_Hull(15);


indbnd = GMMHull.IsInsideHull(Xineq,1.5);
Xineq = Xineq(~indbnd,:);



figure(36)
states=[5,6];
GMMHull.plotGMMpointsHUll(states,Xineq,1.0,'ro')
hold on
plot3(Xineq(:,states(1)),Xineq(:,states(2)),-ones(size(Xineq,1),1),'b+')
hold off
title('GMM HUll')

RIF=RegInterpFitters('DecisionTreeAdaptiveOutRegion');
RIF.fit(dsX.X,dsX.p,Xineq,[],[],GMMHull)


