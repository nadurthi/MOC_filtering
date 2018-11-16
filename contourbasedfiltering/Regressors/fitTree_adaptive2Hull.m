function tree = fitTree_adaptive2Hull(X,p,Pf,Xineq,XtestoutsideHull,GMMhull)
% after you build a tree, generate points in the boxes outside the hull and
% retrain the tree


[N,dim]=size(X);
[Nineq,dimineq]=size(Xineq);

%     factconst = max(p)/10;
%     pnfit = p/factconst;
beq = log(p);

cc=min(beq);
if cc<=0
    bineq = 10*cc*ones(Nineq,1);
else
    bineq = 0.1*(cc+0)*ones(Nineq,1);
end

states=cell(1,dim);
for i=1:dim
    states{i} = num2str(i);
end

if dim==2
    Nmc = 100;
end
if dim == 6
    Nmc = 1000;
end
cnt = 0;

%% augement training data with more points
Xextra=zeros(dim*size(X,1),size(X,2));
bextra = zeros(dim*size(X,1),1);
for i=1:size(X,1)
    Xextra((i-1)*dim+1:i*dim,:) = mvnrnd(X(i,:),0.2^2*eye(dim),dim);
    bextra(i) = beq(i);
end

%% train
Xtrain = [X;Xextra;Xineq];
Btrain = [beq;bextra;bineq];
mean_outof_hullratio = 0;
MAXDIAG=3*sqrt(dim);

lambda =0.7;
boxcutooff = lambda*min(beq)+(1-lambda)*max(beq);

while(1)
    if cnt>3
        break
    end
    cnt=cnt+1;
    disp(cnt)
    tree = fitrtree(Xtrain,Btrain,'MinParentSize',dim+2,'MaxNumSplits',5000,'MinLeafSize',dim+2,...
        'PredictorNames',states,'ResponseName','probs');
    boxes=getTree2Boxes(tree);
    m=0;
    for j=1:size(boxes,1)
        lb = boxes{j,1};
        ub = boxes{j,2};
        y1 = sum(X>repmat(lb,N,1),2)==dim;
        y2 = sum(X<repmat(ub,N,1),2)==dim;
        avgX = mean(beq(y1 & y2));
        
        if avgX < boxcutooff % add points to boxes that have some probability
            continue
        end
        
        maxlen = max(ub-lb);
        diaglen = norm(ub-lb);
        [diaglen,MAXDIAG]
%         if diaglen/MAXDIAG<=0.6
%             disp('box has small diagonal ')
%             continue
%         end
        Xmc=mvurnd(lb,ub,Nmc);
        indbnd = GMMhull.IsInsideHull(Xmc,1.5);
        insidefrac = sum(indbnd)/length(indbnd);
        if insidefrac>0.5
            continue
        end
        maxlen = max(ub-lb);
        
        Xmcout = Xmc(~indbnd,:);
        Xtrain = vertcat(Xtrain,Xmcout);
        Btrain = vertcat(Btrain,bineq(1)*ones(size(Xmcout,1),1));
        
%         m = m + sum(indbnd)/length(indbnd);
    end
    
    disp(['Xtrain has size : ',num2str(size(Xtrain,1))])
    
    mean_outof_hullratio = m;
    
end
end