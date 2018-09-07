function margprobs = get_2Dmarginalized_probs(Xp,keep_n1,keep_n2,X,probs,mX,PX,fullpdf,method)
% Xp is points at which we need probability
% X.probs are the points and weights
% n1,n2 are the dimensions that should remain after marginalization
% fullpdf is corresponding pdf

Xp = Xp(:)';

[xr,xc]=size(X);
if xr==1 || xc==1
    X=X(:)';
end

dim = size(X,2);



% keyboard

keepcols = [keep_n1,keep_n2];
margcols = 1:dim;
margcols(keepcols)=[];
Nmargcols = length(margcols);
Nkeepcols = length(keepcols);
colorder = [margcols,keepcols];

if strcmp(method,'ClusterMC')
    probest = fullpdf.func(X);
    probest=probest/sum(probest);
    
%     m=mquad;
%     Pcov=Pquad;
    
    
    [IDX,C] = kmeans(X, 2);
    remclust = [];
    for i=1:size(C,1)
        [m,pR]=MeanCov(X(IDX==i,:),probest(IDX==i)/sum(probest(IDX==i)));
        if any(eig(pR)<0)
            remclust=[remclust,i];
        end
    end
    C(remclust,:)=[];
    Nclust = size(C,1);
    w = ones(Nclust,1)/Nclust;
    
    Ngh =4;
    Nptscl = Ngh^Nmargcols;
    importpdfeval = zeros(Nptscl,Nclust);
    probinteval = zeros(Nptscl,Nclust);
    gaussclust = cell(1,Nclust);
    for i=1:Nclust
        [m,pR]=MeanCov(X(IDX==i,:),probest(IDX==i)/sum(probest(IDX==i)));
        pR = 1^2*pR;
        mmarg = m(margcols);
        mkeep = m(keepcols);
        
        pRcross = pR(margcols,keepcols);
        pRmarg = pR(margcols,margcols);
        pRkeep = pR(keepcols,keepcols);
        
        [x,wtsq] = GH_points(mmarg,pRmarg,4);
        
%         try
%             x=mvnrnd(mmarg(:)',pRmarg,NMC);
%         catch
%             keyboard
%         end
        
        Y = zeros(Nptscl,dim);
        Y(:,margcols) = x;
        Y(:,keepcols) = repmat(Xp,Nptscl,1);
        
        importpdfeval(:,i) = mvnpdf(x,mmarg(:)',pRmarg);
        probinteval(:,i) = fullpdf.func(Y);
        gaussclust{i} = {m,pR};
    end
%     wtsq = ones(Nptscl,1);
    denompdfval= sum(repmat(w(:)',Nptscl,1).* importpdfeval,2);
    
    
    margprobs = sum( sum( repmat(w(:)',Nptscl,1).*repmat(wtsq(:),1,Nclust).*probinteval,2)./denompdfval );
    
    
    
end


if strcmp(method,'GMM_MC2')
    
    probest = fullpdf.func(X);
    probest=probest/sum(probest);
    
%     m=mquad;
%     Pcov=Pquad;
    
    
    [IDX,C] = kmeans(X, 2);
    remclust = [];
    for i=1:size(C,1)
        [m,pR]=MeanCov(X(IDX==i,:),probest(IDX==i)/sum(probest(IDX==i)));
        if any(eig(pR)<0)
            remclust=[remclust,i];
        end
    end
    C(remclust,:)=[];
    Nclust = size(C,1);
    w = ones(Nclust,1)/Nclust;
    
    NMC=500;
    importpdfeval = zeros(NMC,Nclust);
    probinteval = zeros(NMC,Nclust);
    gaussclust = cell(1,Nclust);
    for i=1:Nclust
        [m,pR]=MeanCov(X(IDX==i,:),probest(IDX==i)/sum(probest(IDX==i)));
        pR = 1^2*pR;
        mmarg = m(margcols);
        mkeep = m(keepcols);
        
        pRcross = pR(margcols,keepcols);
        pRmarg = pR(margcols,margcols);
        pRkeep = pR(keepcols,keepcols);
        
        mmarg_new = mmarg+pRcross*inv(pRkeep)*(Xp(:)-mkeep);
        pRmarg_new = pRmarg-pRcross*inv(pRkeep)*pRcross';
        try
            x=mvnrnd(mmarg_new(:)',pRmarg_new,NMC);
        catch
            keyboard
        end
        
        Y = zeros(NMC,dim);
        Y(:,margcols) = x;
        Y(:,keepcols) = repmat(Xp,NMC,1);
        
        importpdfeval(:,i) = mvnpdf(x,mmarg_new(:)',pRmarg_new);
        probinteval(:,i) = fullpdf.func(Y);
        gaussclust{i} = {m,pR};
    end
    wtsq = ones(NMC,1);
    denompdfval= sum(repmat(w(:)',NMC,1).* importpdfeval,2);
    
    
    margprobs = sum( sum( repmat(w(:)',NMC,1).*repmat(wtsq(:),1,Nclust).*probinteval,2)./denompdfval );
    
    
end




end