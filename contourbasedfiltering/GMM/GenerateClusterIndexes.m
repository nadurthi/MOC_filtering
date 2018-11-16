function idx = GenerateClusterIndexes(X,NgcompMax,method)
[N,dim] = size(X);


previdx=[];
for Ngcomp=3:NgcompMax
    
    if strcmp(method,"kmeans")
        idx = kmeans(X,Ngcomp,'Options',statset('MaxIter',500,'UseParallel',true));
    elseif strcmp(method,"AggClust")
        Z = linkage(X,'complete');
        idx = cluster(Z,'Maxclust',Ngcomp);
    elseif strcmp(method,"gmm")
        GMModel = fitgmdist(X,Ngcomp);
        idx = cluster(GMModel,X) ;
    end
    
    flg=1;
    for i=1:Ngcomp
        xx=X(idx==i,:) ;
        Nc = size(xx,1);
        if Nc <= dim
            disp('BREAK: Nc <= dim')
            flg=0;
            break
        end
        ww = ones(Nc,1)/Nc;
        [mcp,pcp]=MeanCov(xx,ww/sum(ww));
        eigsP = eig(pcp);
        if min(eigsP) <=0 || ~isreal(eigsP) || any(isnan(eigsP)==1) || any(isfinite(eigsP)==0)
            disp('BREAK: min(eigsP) <=0 || isreal(eigsP) || any(isnan(eigsP)==1) || any(isfinite(eigsP)==0)')
            flg=0;
            break
        end
        if cond(pcp)>5000000
            disp('***********\n cond(pcp)>5000000 \n *******************')
            cond(pcp)
            disp(Ngcomp)
            %             flg=0;
            %             break
        end
        sqrt(min(eigsP))
        if sqrt(min(eigsP))<1e-3
            disp('***********\n sqrt(min(eigsP))<1e-3 \n *******************')
%             flg=0;
%             break
        end
    end
    if isempty(previdx)
        previdx = idx;
    end
    if flg==1
        previdx = idx;
    else
        break
    end
    
end

idx = previdx;

end
