function flog = evalKnnPolyFit(X,p,Xineq,Xtest,Pf)
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
    
    Nn=dim+100;
    
    flog=zeros(Ntest,1);
    Xstack = [X;Xineq];
    
    idx = knnsearch(Xstack,Xtest,'K',Nn);
    for i=1:Ntest
        mxentpoly = fitExpPoly_A_Atop_Aineq(Xstack(idx(i,:),:),B(idx(i,:)),Pf,[],[]);
%         mxentpoly = fitleastsquaresPoly(Xstack(idx(i,:),:),B(idx(i,:)),Pf,[],[]);
        
        flog(i) = evaluate_polyND(mxentpoly,Xtest(i,:));
    end
    
end