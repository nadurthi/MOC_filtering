function tree = fitTree(X,p,Pf,Xineq,Xtest)
    [N,dim]=size(X);
    [Nineq,dimineq]=size(Xineq);
    
%     factconst = max(p)/10;
%     pnfit = p/factconst;
    beq = log(p);
    
    cc=min(beq);
    if cc<=0
        bineq = 10*cc*ones(Nineq,1);
    else
        bineq = 0.1*(cc+0)*ones(Nineq,1);
    end

    states=cell(1,dim);
    for i=1:dim
       states{i} = num2str(i); 
    end

    tree = fitrtree([X;Xineq],[beq;bineq],'MinParentSize',dim+2,'MaxNumSplits',5000,'MinLeafSize',dim+2,...
                    'PredictorNames',states,'ResponseName','probs');

end