function [Mdl,B] = fitkdtree(X,p,Xineq)

    [N,dim]=size(X);
    [Nineq,dimineq]=size(Xineq);

    
    beq = log(p);
    
    cc=min(beq);
    if cc<=0
        bineq = 1.5*cc*ones(Nineq,1);
    else
        bineq = 0.1*(cc+0)*ones(Nineq,1);
    end
    B=[beq;bineq];
    
    Nn=dim+1;

    Xtrain  = [X;Xineq];
    Mdl = createns(Xtrain,'BucketSize',Nn);
    