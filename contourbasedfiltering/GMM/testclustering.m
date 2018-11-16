X=mvnrnd([1;1],[3,2.8;2.8,3],100);
plot(X(:,1),X(:,2),'ro')
Z = linkage(X,'complete');

% m = 100
% size of Z is m-1 = 99
% leaf nodes are numbered I = 1 to 100
% each row with index  I is a step that creates a new cluster with cluster id m+I 

NgcompMax=25;
dim=2;
idx = cluster(Z,'Maxclust',NgcompMax);
% [idx,C,sumd,D] = kmeans(X,NgcompMax,'EmptyAction','drop','Options',statset('MaxIter',500,'UseParallel',false));

% [M,P,S]=getmetaCluster(idx,X,NgcompMax);

cnt=0;
while(1)
    cnt=cnt+1;
    cnt
    
    [M,P,S]=getmetaCluster(idx,X,NgcompMax);
    
    B=S(S(:,1)~=-1,:); % take the bad ones
    if isempty(B)
        break
    end
    
    B=B(1,:); % take the first one
    
    badidx = B(1,end);
    RestS = S;
    RestS = RestS(RestS(:,3)~=badidx,:);
    
    badM = M(M(:,3)==badidx,:);
    RestM = M(M(:,3)~=badidx,:);
    
    m=size(RestM,1);
    [~,ind]=min(sum((RestM(:,[1,2])-repmat(badM([1,2]),m,1)).^2,2));
    
    cluttoind = RestM(ind,3);
    
    idx(idx==badidx) = cluttoind;
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
    keyboard
    
end

for i=1:NgcompMax 
    xx=X(idx==i,:) ;

    figure(3)
    plot(X(:,1),X(:,2),'ro')
    hold on
    plot(xx(:,1),xx(:,2),'b+')
    hold off
    
    pause(1)

end


