function margprobs = get_2Dmarginalized_probs2(Xp,n1,n2,X,probs,mX,PX,fullpdf,method)
% Xp is points at which we need probability
% X.probs are the points and weights
% n1,n2 are the dimensions that should remain after marginalization
% fullpdf is corresponding pdf



[xr,xc]=size(Xp);
if xr==1 || xc==1
    Xp=Xp(:)';
end

margprobs = zeros(size(Xp,1),1);


% keyboard

remaincols = [n1,n2];
margcols = 1:size(X,2);
margcols(remaincols)=[];

if strcmp(method,'dummyMC')
%     pp=pdf.func(X);
    pp=probs;
    pp=pp/sum(pp);
%     [m,Pcov]=MeanCov(X(:,margcols),pp);
    m=mX(margcols);
    Pcov=PX(margcols,margcols);
    Pcov = Pcov;
    Nmc=1000;
    

    Xmcfull = zeros(Nmc,size(X,2));
    Xmcfull(:,margcols)=mvnrnd(m,Pcov,Nmc);
    
    for i=1:size(Xp,1)
        Ppartial=get_partial_polyND(P,xfix,fixedinds)
        
        Xmc = Xmcfull;
        Xmc(:,remaincols) = repmat(Xp(i,:),Nmc,1);
        p1 = evaluate_polyND(fullpdf.poly,Xmc);
        S=Xmcfull(:,margcols)-repmat(m(:)',Nmc,1);
        p2=zeros(Nmc,1);
        A=inv(Pcov);
        c1=log(1/sqrt(det(2*pi*Pcov)));
        for j=1:Nmc
            p2(j)=c1-0.5*S(j,:)*A*S(j,:)';
        end
        ppp = p1-p2;
%         ppp(ppp>70)=70;
%         ppp(ppp<-70)=-70;

        
        margprobs(i) = sum(exp(p1-p2).*(1/Nmc));
%         margprobs(i) = sum((fullpdf.func(Xmc)./mvnpdf(Xmcfull(:,margcols),m(:)',Pcov))*1/Nmc);
    end
end





end