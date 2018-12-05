function [M,P,S]=getmetaCluster(idx,X,NgcompMax,mineigvalset)
[N,dim] = size(X);

switch nargin
    case 4
        mineigval=mineigvalset;
    otherwise
        mineigval=1e-3;
end




M=zeros(NgcompMax,dim+1);
P=cell(NgcompMax,1);
S=-1*ones(NgcompMax,2);
for i=1:NgcompMax 
    xx=X(idx==i,:) ;
    Nc = size(xx,1);
    ww = ones(Nc,1)/Nc;

    [m,p]=MeanCov(xx,ww/sum(ww));
    P{i}=p;
    M(i,1:dim)=m;
    M(i,dim+1)=i;
    eigsP = eig(P{i});
    if sqrt(min(eigsP)) <=mineigval || ~isreal(eigsP) || any(isnan(eigsP)==1) || any(isfinite(eigsP)==0)
        disp('BREAK: min(eigsP) <=0 || isreal(eigsP) || any(isnan(eigsP)==1) || any(isfinite(eigsP)==0)')
        S(i,1)=i;
        S(i,2)=real(sqrt(min(eigsP) ));
    end
%     figure(2)
%     plot(X(:,1),X(:,2),'ro')
%     hold on
%     plot(xx(:,1),xx(:,2),'b+')
%     hold off
%     
%     pause(1)
end
S=[S,[1:NgcompMax]'];
