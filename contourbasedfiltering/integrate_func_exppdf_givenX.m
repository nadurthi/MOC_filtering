function intval = integrate_func_exppdf_givenX(F,pdf,X,mquad,Pquad,method)
% X are the points that are tracked
P=pdf.poly;
dim=size(P,2)-1;

%estimate normalizing constant
c=1;
%%  method 1 using MC

if strcmp(method,'dummyMC')
    
    m=mquad;
    Pcov=Pquad;
    
    intval=integrate_func_exppdf(@(x)1,pdf,m,Pcov,method);
    disp(['integration constant : ',num2str(intval)])
    %     keyboard
end

%% method 2 using mixture of gaussians

if strcmp(method,'GMM_MC')
    probest = pdf.func(X);
    probest=probest/sum(probest);
    
%     m=mquad;
%     Pcov=Pquad;
    
    
    [IDX,C] = kmeans(X, 1);
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
    Nptscl = Ngh^dim;
    importpdfeval = zeros(Nptscl,Nclust);
    probinteval = zeros(Nptscl,Nclust);
    gaussclust = cell(1,Nclust);
    for i=1:Nclust
        [m,pR]=MeanCov(X(IDX==i,:),probest(IDX==i)/sum(probest(IDX==i)));
        pR = 0.9^2*pR;
%         x=mvnrnd(m,pR,NMC);
%         Xt=[Xt;x];    
        [x,wtsq] = GH_points(m,pR,4);
        try
            importpdfeval(:,i) = mvnpdf(x,m(:)',pR);
        catch
            keyboard
        end
        probinteval(:,i) = F(x).*pdf.func(x);
        gaussclust{i} = {m,pR};
    end
%     wtsq = ones(NMC,1); 
    denompdfval= sum(repmat(w(:)',Nptscl,1).* importpdfeval,2);
    denompdfval_inv = 1./denompdfval(:);
    
    intval = sum( sum( repmat(w(:)',Nptscl,1).*repmat(wtsq(:),1,Nclust).*probinteval,2)./denompdfval );
    disp(['integration constant : ',num2str(intval)])
    %     keyboard
end




%% method 3 do not use X, but do quadratic approximation to compute mean and covariance, and then do sampling integration






%% now see to append the normalizing constant into the exp pdf







end