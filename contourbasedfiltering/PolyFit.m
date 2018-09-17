
classdef PolyFit < handle
    % fit GMMs and plot the visualizations
    
    properties
        X;
        p;
    end
    methods
        % constructor. save a copy constants
        function obj= PolyFit(X,p)
            % 'TrueState' or 'TransformedState'
            obj.X=X;
            obj.p=p;
        end
        
        function obj = fit(obj)
            factconst = max(pn)/10;
            pnfit = pn/factconst;
            logpnfit = log(pnfit);
            
            
            % [~,ind]=sort(pnfit(:),1,'descend');
            % logppnfit=log(pnfit(ind));
            % AAn=A(ind,:);
            % logppnfit=logppnfit(:);
            
            
            
            % lamdim=length(lam);
            K = (min(logpnfit)-5)*ones(size(Dineq,1),1);
            KK=K;
            DD=Dineq;
            Atop=A(1:15,:);
            logpntop = logpnfit(1:15);
            lenconstr = length(logpnfit);
            
            % %working good
            %     minimize( 10*norm(lam2,1)+50*norm(t,2)+150*norm(t2,2))
            CC=[0.5];
            LAMS=zeros(lamdim,length(CC));
            costs = zeros(1,length(CC));
            for ci = 1:length(CC)
                cvx_begin
                variables t2(15) t(lenconstr) lam2(lamdim)
                minimize( CC(ci)*norm(lam2,1)+50*norm(t,2)+70*norm(t2,2))
                subject to
                %     DD*lam2 <= KK
                A*lam2==logpnfit+t
                Atop*lam2==logpntop+t2
                Tineq*lam2<=max(logpnfit)
                cvx_end
                LAMS(:,ci)=lam2;
                costs(ci) = norm(t,2);
            end
            [~,bind] = min(costs);
            lam2 = LAMS(:,bind);
            lam2(1) = lam2(1)+log(factconst);
            lamsol = lam2;
            
            mxentpoly_norm=zeros(1,dim+1);
            for i=1:length(lamsol)
                mxentpoly_norm=add_sub_polyND(mxentpoly_norm, scalar_multiply_polyND(lamsol(i),Pf{i}),'add');
            end
            mxentpoly_norm=simplify_polyND(mxentpoly_norm);
            
            pdfnorm.func=@(x)exp(evaluate_polyND(mxentpoly_norm,x));
            pdfnorm.polyeval=@(x)evaluate_polyND(mxentpoly_norm,x);
            pdfnorm.poly=mxentpoly_norm;
            pdfnorm.type = 'true-0I-hypercube-11';
            
            
            
        end
        
    end
end