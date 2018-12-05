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
if dim == 4
    Nmc = 300;
end

if dim == 6
    Nmc = 1000;
end


%% augement training data with more points
NextraMC=5;
Nlimit=50;
Xextra=zeros(NextraMC*Nlimit,dim);
bextra = zeros(NextraMC*Nlimit,1);
for i=1:Nlimit
    Xextra((i-1)*NextraMC+1:i*NextraMC,:) = mvnrnd(X(i,:),0.05^2*eye(dim),NextraMC);
    bextra(i) = beq(i);
end

%% train
Xtrain = [X;Xextra;Xineq];
Btrain = [beq;bextra;bineq];
mean_outof_hullratio = 0;
MAXDIAG=3*sqrt(dim);

lambda =0.7;
boxcutooff = lambda*min(beq)+(1-lambda)*max(beq);

% keyboard

cnt = 0;
while(1)
    if cnt>=2
        break
    end
    cnt=cnt+1;
    disp(cnt)
    tree = fitrtree(Xtrain,Btrain,'MinParentSize',dim+2,'MaxNumSplits',5000,'MinLeafSize',dim+2,...
        'PredictorNames',states,'ResponseName','probs');
    tic
    boxes=getTree2Boxes(tree);
    toc
    
    m=0;
    XX=zeros(100000,dim);
    k=1;
    for j=1:size(boxes,1)
        [j,size(boxes,1)]
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
%         [diaglen,MAXDIAG]
%         if diaglen/MAXDIAG<=0.6
%             disp('box has small diagonal ')
%             continue
%         end
        Xmc=mvurnd(lb,ub,Nmc);
        indbnd = GMMhull.IsInsideHull(Xmc,1.15);
        insidefrac = sum(indbnd)/length(indbnd);
        if insidefrac>0.5
            continue
        end
        maxlen = max(ub-lb);
        
        Xmcout = Xmc(~indbnd,:);
        Nmcout  = size(Xmcout,1);
        XX(k:k+Nmcout-1,:)=Xmcout;
        k = k+Nmcout;
%         m = m + sum(indbnd)/length(indbnd);
    end

        XXX = XX(1:k-1,:);
        Xtrain = vertcat(Xtrain,XXX);
        Btrain = vertcat(Btrain,bineq(1)*ones(size(XXX,1),1));

    disp(['Xtrain has size : ',num2str(size(Xtrain,1))])
    
%     mean_outof_hullratio = m;
    
end
end