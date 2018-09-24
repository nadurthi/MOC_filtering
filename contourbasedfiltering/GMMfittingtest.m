dsX = DataSet(X,probs,'TrueState');
[N,dim] =size(X);

GMModel = fitgmdist(X,50);
idx = cluster(GMModel,X) ;
        
% dsX.AddHyperCubeTrasform(-1*ones(dim,1),1*ones(dim,1));

dsX = DataSet(X,probs,'TrueState');
[N,dim] =size(X);

mX=mX(:);

dsX.AddMeanCov_to_OI_Trasform(mX,PX);
dsX.AddHyperCubeTrasform(-1*ones(dim,1),1*ones(dim,1));

figure(33)
dsX.PlotPointsProbs3D([1,2],'ro');
title(['true points and MC: time step = ',num2str(4)])
hold off

[idx,C]=kmeans(dsX.X,500);
figure
for i=1:100
    clf
    plot(C(:,1),C(:,2),'k*')
    hold on
   plot(dsX.X(idx==i,1),dsX.X(idx==i,2),'ro')
   axis([-1,1,-1,1])
   pause(0.2)
end

GMMfitter = GMMFitDataSet(dsX.X,dsX.p);
GMMfitter.SetGMM_Hull(200);
figure(36)
GMMfitter.plotGMMpointsHUll([1,2],2,'ro')
hold on
plot3(Xineq(:,1),Xineq(:,2),-ones(size(Xineq,1),1),'b+')
hold off
title('GMM HUll')

% GMM = GMMfitter.FitGMM_1comp();
GMM = GMMfitter.FitGMM_kmeans_optimwt(3);

% figure(36)
% GMMfitter.plotGMMpointsHUll([1,2],2,'ro')

figure(34)
GMMfitter.plotGMMpoints([1,2],'ro')
title('GMM fit points')
hold off

figure(35)
GMMfitter.plotGMMSurf([1,2],'ro')
title('GMM fit surf')
hold off
Xineq=[];
% [Xbnd1,~] = GLgn_pts(-1.5*ones(1,dim),1.5*ones(1,dim),5);

SS = mvurnd(-2*ones(dim,1),2*ones(dim,1),10000);

% Xbnd2=2.5*gen_uniform_grid(3,dim);
% Xbnd2=5*gen_uniform_grid(5,dim);
Xineq = [SS];

GMMfitter.SetGMM_Hull();
indbnd = GMMfitter.IsInsideHull(Xineq,2);
Xineq = Xineq(~indbnd,:);
figure(36)
GMMfitter.plotGMMpointsHUll([1,2],2,'ro')
hold on
plot3(Xineq(:,1),Xineq(:,2),-ones(size(Xineq,1),1),'b+')
hold off
title('GMM HUll')
Nineq = size(Xineq,1);

x = [dsX.X(:,1);Xineq(:,1)];
y = [dsX.X(:,2);Xineq(:,1)];
z = [dsX.X(:,3);Xineq(:,1)];
vx = [dsX.X(:,4);Xineq(:,1)];
vy = [dsX.X(:,5);Xineq(:,1)];
vz = [dsX.X(:,6);Xineq(:,1)];

Tbl = table(x,y,z,vx,vy,vy);

tree = fitrtree(Tbl,[log(dsX.p);min(log(dsX.p))*ones(Nineq,1)],'MinParentSize',5,'MinLeafSize',5);

SS = mvurnd(-1*ones(dim,1),1*ones(dim,1),10000);

py=predict(tree,SS);

figure
plot3(dsX.X(:,1),dsX.X(:,2),log(dsX.p),'ro')
hold on
plot3(SS(:,1),SS(:,2),py,'b+')

