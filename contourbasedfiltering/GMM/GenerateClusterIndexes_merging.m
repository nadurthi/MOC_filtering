function idx = GenerateClusterIndexes_merging(X,NgcompMax,method,mineigvalset)
[N,dim] = size(X);

switch nargin
    case 4
        mineigval=mineigvalset;
    otherwise
        mineigval=1e-3;
end


Z = linkage(X,'complete');
idx = cluster(Z,'Maxclust',NgcompMax);

cnt=0;
while(1)
    cnt=cnt+1;
    cnt
    
    [M,P,S]=getmetaCluster(idx,X,NgcompMax,mineigval);
    
    B=S(S(:,1)~=-1,:); % take the bad ones
    if isempty(B)
        break
    end
    
    B=B(1,:); % take the first one
    
    badidx = B(1,end);
    RestS = S;
    RestS = RestS(RestS(:,3)~=badidx,:);
    
    badM = M(M(:,end)==badidx,:);
    RestM = M(M(:,end)~=badidx,:);
    
    m=size(RestM,1);
    [~,ind]=min(sum((RestM(:,1:dim)-repmat(badM(1:dim),m,1)).^2,2));
    
    cluttoind = RestM(ind,end);
    if isempty(cluttoind)
        break
    end
    try
        idx(idx==badidx) = cluttoind;
    catch
        keyboard
    end
    if cluttoind>badidx
        for i=badidx+1:length(unique(idx))+1
           idx(idx==i)=idx(idx==i)-1; 
        end
    else
        for i=badidx+1:length(unique(idx))+1
           idx(idx==i)=idx(idx==i)-1; 
        end
    end
    % make the idex cosecutinve
%     if length(unique(idx))~=max(idx)
%         h=unique(idx);
%         
%     end
    [length(unique(idx)),max(idx)]
    NgcompMax = max(idx);
%     keyboard
    
end

% idx = previdx;

end