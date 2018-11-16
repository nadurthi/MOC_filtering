function idx = GenerateClusterIndexes_merging2222(X,NgcompMax,method)
[N,dim] = size(X);

if strcmp(method,"kmeans")
    idx = kmeans(X,NgcompMax,'Options',statset('MaxIter',500,'UseParallel',true));
elseif strcmp(method,"AggClust")
    Z = linkage(X,'complete');
    idx = cluster(Z,'Maxclust',NgcompMax);
elseif strcmp(method,"gmm")
    GMModel = fitgmdist(X,NgcompMax);
    idx = cluster(GMModel,X) ;
end

while(1)
    M=zeros(NgcompMax,dim+1);
    P=cell(NgcompMax,1);
    S=ones(NgcompMax,1);
    for i=1:NgcompMax 
        xx=X(idx==i,:) ;
        Nc = size(xx,1);
        ww = ones(Nc,1)/Nc;

        [m,p]=MeanCov(xx,ww/sum(ww));
        P{i}=p;
        M(i,1:dim)=m;
        M(i,dim+1)=i;
        eigsP = eig(P{i});
        if sqrt(min(eigsP)) <=1e-3 || ~isreal(eigsP) || any(isnan(eigsP)==1) || any(isfinite(eigsP)==0)
            disp('BREAK: min(eigsP) <=0 || isreal(eigsP) || any(isnan(eigsP)==1) || any(isfinite(eigsP)==0)')
            S(i)=0;
        end
    end
    
    if all(S==1)
       break
    end

    ind = S==0;
    badM = M(ind,:);
    Mrest = M(M(:,dim+1)~=badM(1,dim+1),:);
    gg = knnsearch(Mrest,badM(1,:),'K',1);
    
    minidx = min([gg,badM(1,end)]);
    idx(idx==gg)=minidx;
    idx(idx==badM(1,end))=minidx;

    NgcompMax = max(idx);
    [gg,badM(1,end),NgcompMax]
end


% idx = previdx;

end