function [priorpdfk1,postpdfk1,Xnptk1_prior,probntk1_prior,Xnptk1_post,probntk1_post]=moc_pf_6D(pdfk,tk,time,zk1,model,optns)

% if zk is Nan do not do any measurement update
% pdfk.normeval, pdfk.tree, [pdfk.Xrep,pdfk.prep],  ... optional [pdfk.X, pdfk.p]
% tk,tk1 are the steps
% tk1 = tk+1
% first use UKF to propagate from the backk  to k1, using reppoints Xrep
% compute mean and cov, take 5sigma, go back and get more points
%%
% [pdfk.Xrep,pdfk.prep] are real-space points in
[N,dim] =size(pdfk.Xrep);

[Xrepk1,prepk1]=propagate_character(pdfk.Xrep,pdfk.prep,time.dt,time.Tvec(tk),model);

[mXk1,PXk1]=MeanCov(Xrepk1,prepk1/sum(prepk1));


dsX = DataSet(Xrepk1,prepk1,'TrueState');


% mXk1=Xk1(b,:)';
dsX.AddMeanCov_to_OI_Trasform(mXk1,3^2*PXk1);
% dsX.AddHyperCubeTrasform(-0.9*ones(dim,1),0.9*ones(dim,1));

plot3(dsX.X(:,1),dsX.X(:,2),dsX.p,'ro')

% keyboard

% indd = dsX.p<1e-100;
% dsX.p(indd) = 1e-100;

dsX.SortByProb('descend');


%% plottinmg
% for i=0:0
%     figure(80+i)
%     dsX.PlotPointsProbs3D([i+1,i+2],'ro');
% end


%% Form the hyepercube and go get back the fucntion values for prior pdf


LB=-1.4*ones(dim,1);
UB=1.4*ones(dim,1);

if optns.genpriorpdf==1
    Xnptk1=[];
    probntk1=[];
    cnt=0;
    while 1
        % first get the test points
        if isempty(Xnptk1)==1
%             XXnptk1 = mvurnd(1*LB,1*UB,dim*100);
            XXnptk1=[dsX.X;mvurnd_expanded(1*LB,1*UB,dim*1000,0)];
%             probntk1=dsX.p
        else
            % tree based adaptive sampling
            
            nb=size(boxes,1);
            %         Pbox=zeros(nb,1);
            %         for i=1:nb
            %            Pbox(i)=boxes{i,4};
            %         end
            %         Pbox=Pbox/sum(Pbox);
            

%             nids = discretesample(errbox, dim*500);
            nids = pmf_sampler(Mbox,dim*500);
%             nids = [setdiff(1:nb,nids(:)'),nids(:)'];

            XXnptk1=[];
            for i=1:nb
                ns = sum(nids==i);
                lb = boxes{i,1};
                ub = boxes{i,2};
                XXnptk1=[XXnptk1;mvurnd_expanded(lb,ub,2*ns,0)];
            end
            
            nids=pmf_sampler(Assratio,dim*50);
            for i=1:nb
                ns = sum(nids==i);
                lb = boxes{i,1};
                ub = boxes{i,2};
                XXnptk1=[XXnptk1;mvurnd_expanded(lb,ub,2*ns,0)];
            end
        end
        Nm=size(XXnptk1,1);
        [XXptk1,~] = dsX.ApplyAffineTransform_Final2Original(XXnptk1,ones(Nm,1));
        [XXptk,~] = propagate_character_back(XXptk1,ones(Nm,1),time.dt,time.Tvec(tk+1),model);
        
        XXnptk=pdfk.transForms.trueX2normX(XXptk);
        Pprobntk = pdfk.normeval(XXnptk);
        
        Pprobtk=pdfk.transForms.normprob2trueprob(Pprobntk);
        % for now we take trace of dynamics is 0
        Pprobtk1 = Pprobtk;
        
        Pprobntk1 = dsX.ApplyAffineTransformProb_Original2Final(Pprobtk1);
        
        figure(33)
        plot3(XXnptk1(:,optns.plotstates(1)),XXnptk1(:,optns.plotstates(2)),Pprobntk1,'bo')
        
        Xnptk1=[Xnptk1;XXnptk1];
        probntk1=[probntk1;Pprobntk1];
        
        Nptk1 = size(Xnptk1,1);
        
        % [Xnptk1,probntk1];
        Ntree=1;
        EE=cell(Ntree,3);
        disp('EE')
        tic
        for jj=1:Ntree
            cv = cvpartition(size(Xnptk1,1),'HoldOut',0.3);
            idx = cv.test;

            Xtrain=Xnptk1(~idx,:);
            ytrain=probntk1(~idx);
            Xtest=Xnptk1(idx,:);
            ytest=probntk1(idx);

            % Now build the tree
            tree = fitrtree(Xtrain,ytrain,'MinParentSize',6,'MaxNumSplits',5000,'MinLeafSize',6);
            ypred=tree.predict(Xtest);
            errtest=mean(100*abs(ypred-ytest)./(ytest+1));
            EE{jj,1}=tree;
            EE{jj,2}=errtest;
            
            boxes=getTree2Boxes(tree,LB,UB);
            AA=zeros(size(boxes,1),1);
            for r=1:size(boxes,1)
                lb = boxes{r,1};
                ub = boxes{r,2};
                d=ub-lb;
                AA(r) = abs(max(d)/min(d));
                if isnan(AA(r)) || isinf(AA(r))
                    AA(r)=1e3;
                end
                if AA(r)==0
                    AA(r)=1;
                end
                if AA(r)>1e10
                    AA(r)=1e9;
                end
            end
            EE{jj,3}=max(AA);
        end
        toc
        EE=sortrows(EE,2);
        EE=EE(1:1,:);
        EE=sortrows(EE,3);
        tree=EE{1,1};
        errtest=EE{1,2};
        disp('box')

        
        boxes=getTree2Boxes(tree,LB,UB);
        boxes = sortrows(boxes,4,'descend');
        
        disp('plot')
        tic
        figure(71)
        plot3(Xnptk1(:,optns.plotstates(1)),Xnptk1(:,optns.plotstates(2)),probntk1,'bo')
        hold on
        plot3(Xnptk1(:,optns.plotstates(1)),Xnptk1(:,optns.plotstates(2)),tree.predict(Xnptk1),'r+')
        plot3Dtreepatches_6D(boxes(1:20,:),optns.plotstates)
        hold off
        toc

        
        errbox=zeros(size(boxes,1),1);
        Mbox=zeros(size(boxes,1),1);
        Assratio=zeros(size(boxes,1),1);
        disp('Mbox')
        tic

            
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
            y1 = sum(Xnptk1>repmat(lb(:)',Nptk1,1),2)==dim;
            y2 = sum(Xnptk1<repmat(ub(:)',Nptk1,1),2)==dim;
            a = probntk1(y1 & y2);
            m = boxes{j,4};
            if isempty(a)==1
                errbox(j)=0;
            else
                errbox(j) = mean(abs(a-repmat(m,length(a),1))/(m+1));
            end
            Mbox(j) = m;
        end
        if any(isnan(errbox))
            keyboard
        end

        toc
        
        disp('s%')
        tic
        s=prctile(Mbox,65);
        inds = Mbox>s;
        errbox(inds)=errbox(inds)*5;
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
        Assratio(inds)=Assratio(inds)*5;
        if sum(Assratio)~=0
            Assratio = Assratio/sum(Assratio);
        end
        toc
        

        
%         keyboard
                
        
        if (errtest<8 && cnt>5) || cnt>6
            break
        end
        cnt=cnt+1;
        [cnt,errtest,Nptk1]
        
    end
    
    Xnptk1_prior = Xnptk1;
    probntk1_prior = probntk1;
    
    priorpdfk1.normeval=@(x)treepredict_withbnds(x,tree,LB,UB);
    priorpdfk1.tree=tree;
    priorpdfk1.boxes=boxes;
    priorpdfk1.dim =dim;
    priorpdfk1.transForms = dsX.GetTrasnformers();
    priorpdfk1.info = 'prior-pdf-tree';
    priorpdfk1.pdftype = 'tree';
    priorpdfk1.LB = LB;
    priorpdfk1.UB = UB;
    
    nb=size(boxes,1);
    nids=pmf_sampler(Mbox/sum(Mbox),dim*100);
    XX=[];
    for i=1:nb
        ns = sum(nids==i);
        lb = boxes{i,1};
        ub = boxes{i,2};
        XX=[XX;mvurnd_expanded(lb,ub,10*ns,0)];
    end
    
    for i=1:nb
        lb = boxes{i,1};
        ub = boxes{i,2};
        XX=[XX;mvurnd_expanded(lb,ub,5,0)];
    end
    
    pp = priorpdfk1.normeval(XX);
    [priorpdfk1.Xrep,priorpdfk1.prep] = dsX.ApplyAffineTransform_Final2Original(XX,pp);

    
else
    probntk1_prior = NaN;
    Xnptk1_prior = NaN;
    priorpdfk1 = NaN;
end
disp('-------------------')
tic
figure(71)
plot3(Xnptk1(:,1),Xnptk1(:,2),probntk1,'bo')
hold on
plot3(Xnptk1(:,1),Xnptk1(:,2),tree.predict(Xnptk1),'r+')
% plot3Dtreepatches_6D(boxes,[1,2])
plot2Dboxes_modf6D(boxes,[1,2])
hold off

figure(72)
plot3(Xnptk1(:,2),Xnptk1(:,3),probntk1,'bo')
hold on
plot3(Xnptk1(:,2),Xnptk1(:,3),tree.predict(Xnptk1),'r+')
plot2Dboxes_modf6D(boxes,[2,3])
hold off

figure(73)
plot3(Xnptk1(:,3),Xnptk1(:,4),probntk1,'bo')
hold on
plot3(Xnptk1(:,3),Xnptk1(:,4),tree.predict(Xnptk1),'r+')
plot2Dboxes_modf6D(boxes,[3,4])
hold off

figure(74)
plot3(Xnptk1(:,4),Xnptk1(:,5),probntk1,'bo')
hold on
plot3(Xnptk1(:,4),Xnptk1(:,5),tree.predict(Xnptk1),'r+')
plot2Dboxes_modf6D(boxes,[4,5])
hold off

figure(75)
plot3(Xnptk1(:,5),Xnptk1(:,6),probntk1,'bo')
hold on
plot3(Xnptk1(:,5),Xnptk1(:,6),tree.predict(Xnptk1),'r+')
plot2Dboxes_modf6D(boxes,[5,6])
hold off

toc

% keyboard
%% Now update posterior pdf directly from pdfk

if optns.genpostpdf==1
    Xnptk1=[];
    probntk1=[];
    cnt=0;
    while 1
        % first get the test points
        % first get the test points
        if isempty(Xnptk1)==1
%             XXnptk1 = mvurnd(1*LB,1*UB,dim*100);
            XXnptk1=[dsX.X;mvurnd_expanded(1*LB,1*UB,dim*50,0)];
%             probntk1=dsX.p
        else
            % tree based adaptive sampling
            
            nb=size(boxes,1);
            %         Pbox=zeros(nb,1);
            %         for i=1:nb
            %            Pbox(i)=boxes{i,4};
            %         end
            %         Pbox=Pbox/sum(Pbox);
            
            try
            nids = discretesample(errbox, dim*50);
            nids = [setdiff(1:nb,nids(:)'),nids(:)'];
            catch
                keyboard
            end
            XXnptk1=[];
            for i=1:nb
                ns = sum(nids==i);
                lb = boxes{i,1};
                ub = boxes{i,2};
                XXnptk1=[XXnptk1;mvurnd_expanded(lb,ub,5*ns,0)];
            end
            
            nids=pmf_sampler(Assratio,dim*10);
            for i=1:nb
                ns = sum(nids==i);
                lb = boxes{i,1};
                ub = boxes{i,2};
                XXnptk1=[XXnptk1;mvurnd_expanded(lb,ub,10*ns,0)];
            end
        end
        Nm=size(XXnptk1,1);
        [XXptk1,~] = dsX.ApplyAffineTransform_Final2Original(XXnptk1,ones(Nm,1));
        [XXptk,~]=propagate_character_back(XXptk1,ones(Nm,1),time.dt,time.Tvec(tk+1),model);
        
        XXnptk=pdfk.transForms.trueX2normX(XXptk);
        Pprobntk = pdfk.normeval(XXnptk);
        
        Pprobtk=pdfk.transForms.normprob2trueprob(Pprobntk);
        % for now we take trace of dynamics is 0
        Pprobtk1 = Pprobtk;
        
        % do meas update with zk
        logpz = log(1);

        Pprobtk1(Pprobtk1<1e-70)=1e-70;
        
        logprobs=log(Pprobtk1);
        logprobsXpost = zeros(size(Pprobtk1));

        for i=1:size(XXptk1,1)
            logprobsXpost(i) = log(1/sqrt(det(2*pi*model.R)))-0.5*(zk1(:)-model.h(XXptk1(i,:)'))'*inv(model.R)*(zk1(:)-model.h(XXptk1(i,:)'))+logprobs(i)-logpz;
        end

        Pprobtk1post = exp(logprobsXpost);

        Pprobntk1 = dsX.ApplyAffineTransformProb_Original2Final(Pprobtk1post);
        
        Xnptk1=[Xnptk1;XXnptk1];
        probntk1=[probntk1;Pprobntk1];
        
        Nptk1 = size(Xnptk1,1);
        
        % [Xnptk1,probntk1];
        EE=cell(25,2);
        disp('EE')
        tic
        for jj=1:25
            cv = cvpartition(size(Xnptk1,1),'HoldOut',0.3);
            idx = cv.test;

            Xtrain=Xnptk1(~idx,:);
            ytrain=probntk1(~idx);
            Xtest=Xnptk1(idx,:);
            ytest=probntk1(idx);

            % Now build the tree
            tree = fitrtree(Xtrain,ytrain,'MinParentSize',dim*2,'MaxNumSplits',5000,'MinLeafSize',dim*2);
            ypred=tree.predict(Xtest);
            errtest=max(100*abs(ypred-ytest)./(ytest+1));
            EE{jj,1}=tree;
            EE{jj,2}=errtest;
        end
        toc
        EE=sortrows(EE,2);
        tree=EE{1,1};
        errtest=EE{1,2};
        disp('box')
        tic
        boxes=getTree2Boxes(tree,LB,UB);
        toc
        boxes = sortrows(boxes,4,'descend');
        
%         figure(82)
%         plot3(Xnptk1(:,1),Xnptk1(:,2),probntk1,'bo')
%         hold on
%         plot3Dtreepatches(boxes)
%         hold off
        
%         
        
        errbox=zeros(size(boxes,1),1);
        Mbox=zeros(size(boxes,1),1);
        Assratio=zeros(size(boxes,1),1);
        disp('Mbox')
        tic

            
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
            y1 = sum(Xnptk1>repmat(lb(:)',Nptk1,1),2)==dim;
            y2 = sum(Xnptk1<repmat(ub(:)',Nptk1,1),2)==dim;
            a = probntk1(y1 & y2);
            m = boxes{j,4};
            if isempty(a)==1
                errbox(j)=0;
            else
                errbox(j) = mean(abs(a-repmat(m,length(a),1))/(m+1));
            end
            Mbox(j) = m;
        end
        if any(isnan(errbox))
            keyboard
        end

        toc
        
        disp('s%')
        tic
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
        toc
        

        
%         keyboard
                
        
        if (errtest<25 && cnt>5) || cnt>25
            break
        end
        cnt=cnt+1;
        [cnt,Nptk1]
        
    end
    % get normalization constant
    V=0;
    for j=1:size(boxes,1)
        lb = boxes{j,1};
        ub = boxes{j,2};
        m = boxes{j,4};
        V=V+prod(abs(ub-lb))*m;
    end
    for j=1:size(boxes,1)
        boxes{j,4}=boxes{j,4}/V;
    end
    
    probntk1_post = probntk1;
    Xnptk1_post = Xnptk1;
    
    postpdfk1.normeval=@(x)(1/V)*treepredict_withbnds(x,tree,LB,UB);
    postpdfk1.volumeconstant = V;
    postpdfk1.tree=tree;
    postpdfk1.boxes=boxes;
    postpdfk1.dim =dim;
    postpdfk1.transForms = dsX.GetTrasnformers();
    postpdfk1.info = 'post-pdf-tree';
    postpdfk1.pdftype = 'tree';
    postpdfk1.LB = LB;
    postpdfk1.UB = UB;
    
    % generate representitative points
    nb=size(boxes,1);
    nids=pmf_sampler(Mbox/sum(Mbox),dim*75);
    nids=[nids(:)',setdiff(1:nb,nids(:)')];
    XX=[];
    for i=1:nb
        ns = sum(nids==i);
        lb = boxes{i,1};
        ub = boxes{i,2};
        XX=[XX;mvurnd_expanded(lb,ub,10*ns,0)];
    end
    pp = postpdfk1.normeval(XX);
    [postpdfk1.Xrep,postpdfk1.prep] = dsX.ApplyAffineTransform_Final2Original(XX,pp);

    figure(82)
    plot3(Xnptk1(:,optns.plotstates(1)),Xnptk1(:,optns.plotstates(2)),probntk1,'bo')
    hold on
    plot3Dtreepatches_6D(boxes,optns.plotstates)
    hold off
    toc

else
    probntk1_post = NaN;
    Xnptk1_post = NaN;
    postpdfk1=NaN;
end

