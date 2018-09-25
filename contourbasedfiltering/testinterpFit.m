dsX = DataSet(X,probs,'TrueState');
[N,dim] =size(X);


dsX.AddMeanCov_to_OI_Trasform(mX,PX);
dsX.AddHyperCubeTrasform(-1*ones(dim,1),1*ones(dim,1));

indd = dsX.p<1e-70;
dsX.p(indd) = 1e-70;

dsX.SortByProb('descend');

logpn = dsX.getlogProb();

clc
close all
figure(33)
dsX.PlotPointsProbs3D([1,2],'ro');
figure(34)
dsX.PlotPointsProbs3D([2,3],'ro');
figure(35)
dsX.PlotPointsProbs3D([3,4],'ro');
figure(36)
dsX.PlotPointsProbs3D([4,5],'ro');
figure(37)
dsX.PlotPointsProbs3D([5,6],'ro');

%%
clc
LB=-1.5*ones(dim,1);
UB=1.5*ones(dim,1);

LBtest=-1.5*ones(dim,1);
UBtest=1.5*ones(dim,1);

% XtestingMC=1.5*gen_uniform_grid(5,6);
XtestingMC = mvurnd(LB,UB,500000);
% XtestingMC = dsX.X;

 idx = knnsearch(dsX.X,XtestingMC,'K',5,'Distance','cityblock');
%   idxtest = knnsearch(dsX.X,xtest,'K',5,'Distance','cityblock');
 
  
 pmc = mean(dsX.p(idx),2);
 
%  pmc = exp(logpmc);

 states = [1,2]
figure(33)
plot3(dsX.X(:,states(1)),dsX.X(:,states(2)),dsX.p,'ro')
hold on
% plot3(dsX.X(idxtest,states(1)),dsX.X(idxtest,states(2)),dsX.p(idxtest),'ks','linewidth',2)
plot3(XtestingMC(:,states(1)),XtestingMC(:,states(2)),pmc,'b+')
% Vmc = interpn(dsX.X(:,1),dsX.X(:,2),dsX.X(:,3),dsX.X(:,4),dsX.X(:,6),dsX.X(:,6),logpn,XtestingMC(:,1),XtestingMC(:,2),XtestingMC(:,3),XtestingMC(:,4),XtestingMC(:,5),XtestingMC(:,6));
hold off

