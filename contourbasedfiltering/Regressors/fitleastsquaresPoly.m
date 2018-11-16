function mxentpoly = fitleastsquaresPoly(X,p,Pf,Xineq,Xtest)
    [N,dim] = size(X);
    lamdim = length(Pf);
    Aeq=zeros(N,lamdim);
    

    for ib=1:length(Pf)
        Aeq(:,ib)=evaluate_PolyforLargetSetX(Pf{ib},X);
    end



    factconst = max(p)/10;
    pnfit = p/factconst;
    beq = log(pnfit);

    lam = ( (Aeq'*Aeq)\Aeq')*beq;

    lam(1) = lam(1)+log(factconst);
    lamsol = lam;

    mxentpoly=zeros(1,dim+1);
    for i=1:lamdim
        mxentpoly=add_sub_polyND(mxentpoly, scalar_multiply_polyND(lamsol(i),Pf{i}),'add');
    end
    mxentpoly=simplify_polyND(mxentpoly);


end