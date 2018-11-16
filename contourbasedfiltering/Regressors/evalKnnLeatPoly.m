function f = evalKnnLeatPoly(X,p,Xineq,Xtest)
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
    
    Nn=dim+5;
    
    f=zeros(Ntest,1);
    Xstack = [X;Xineq];
    
    idx = knnsearch(Xstack,Xtest,'K',Nn);
    for i=1:Ntest
        [mp,mpfactor] = mpapi(Xstack(idx(i,:),:)',B(idx(i,:))');
        f(i) = mpval(mp,Xtest(i,:)');
    end
    
end