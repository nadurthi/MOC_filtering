function f = evalKnnMean(X,p,Xineq,Xtest)
    [N,dim]=size(X);
    [Nineq,dimineq]=size(Xineq);
    [Ntest,dimtest]=size(Xtest);
    
    beq = log(p);
    
    cc=min(beq);
    if cc<=0
        bineq = 1.5*cc*ones(Nineq,1);
    else
        bineq = 0.1*(cc+0)*ones(Nineq,1);
    end
    B=[beq;bineq];
    
    Nn=dim+1;
%     disp(['in evalknnmean'])
    f=zeros(Ntest,1);
    Xtrain  = [X;Xineq];
    idx = knnsearch(Xtrain,Xtest,'K',Nn);
%     Xtrain(idx)
%     keyboard
    for i=1:Ntest
        dists = sqrt(sum((Xtrain(idx(i,:),:)-repmat(Xtest(i,:),Nn,1)).^2,2));
        wts = 1./dists;
        wts = wts/sum(wts);
        f(i) = sum(wts.*B(idx(i,:)));
    end

           
end