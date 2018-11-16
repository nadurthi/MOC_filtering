function [GMM,boxes] = fitGMM2tree(X,p,tree)
    [N,dim]=size(X);
    boxes=getTree2Boxes(tree);
    
    Nb = size(boxes,1);
    
    GMM.w = ones(Nb,1)/Nb;
    GMM.mx = cell(Nb,1);
    GMM.Px = cell(Nb,1);
    GMM.Ngcomp = Nb;

    for i=1:Nb
        lb = boxes{i,1};
        ub = boxes{i,2};
        y1 = sum( (X - repmat(lb,N,1))>=0,2);
        y2 = sum( (repmat(ub,N,1)-X)>=0,2);
        
        Xd = X((y1==dim) & (y2==dim),:);
        pd = p((y1==dim) & (y2==dim));

        
%         # NEED TO COMPUTE THE     
        if isempty(Xd)
            Xd=mvurnd(lb,ub,100);
            pd=0.2*min(p)*ones(100,1);
        end
        
        [mm,PP] = MeanCov(Xd,pd/sum(pd));
        [U,D]=eig(PP);
        if any(diag(D)<=1e-4)
            dd=diag(D);
            dd(dd<1e-4)=1e-4;
            D=diag(dd);
            PP=U*D*U'; 
        end

        
        ww=GMM.w(i);
%         [c,err,flag] = fmincon(@(c)norm(pd-ww*mvnpdf(Xd,mm(:)',c*PP)),1,[],[],[],[],1,1e10);
%         if flag<0
%             keyboard
%         end
        c=1;
        
        GMM.mx{i} = mm;
        GMM.Px{i} = c*PP;
    end
%     keyboard
    % now optimize weights
    A=zeros(N,GMM.Ngcomp);
    b=zeros(N,1);
    
    b=p(1:N);
    for i=1:GMM.Ngcomp
       A(:,i) =  mvnpdf(X(1:N,:),GMM.mx{i}',GMM.Px{i});
    end
    ncomp = GMM.Ngcomp;
    cvx_begin
        variables w(ncomp) tt(N)
        minimize( norm(tt) )
        subject to
        A*w==b+tt;
        w>=0;
        w<=1;
        sum(w)<=1.2;
        sum(w)>=0.8;
    cvx_end
    
    GMM.w = w;
    
end