
classdef PolyFit < handle
    % fit GMMs and plot the visualizations
    
    properties
        X;
        p;
        poly;
        PfBasis;
        mxentpoly;
    end
    methods
        % constructor. save a copy constants
        function obj= PolyFit(X,p)
            obj.X=X;
            obj.p=p;
        end

        function mxentpoly = fitExpPoly_A_Atop_Aineq(obj,Pf,Xineq)
            [N,dim] = size(obj.X);
            lamdim = length(Pf);
            Aeq=zeros(N,lamdim);
            for r=1:1:N
                Aeq(r,:) = evaluate_MatrixOfPolys(Pf,obj.X(r,:));
            end

            factconst = max(obj.p)/10;
            pnfit = obj.p/factconst;
            beq = log(pnfit);
            
            Ntop = 15;
            AeqTop = Aeq(1:Ntop,:);
            beqTop = beq(1:Ntop);
            
            % lamdim=length(lam);
            Aineq=zeros(size(Xineq,1),lamdim);
            for r=1:1:size(Xineq,1)
                Aineq(r,:) = evaluate_MatrixOfPolys(Pf,Xineq(r,:));
            end
            bineq = (min(beq)-0 )*ones(size(Aineq,1),1);
            
            
            
            % %working good
            %     minimize( 10*norm(lam2,1)+50*norm(t,2)+150*norm(t2,2))
            CC=[0.005];
            LAMS=zeros(lamdim,length(CC));
            costs = zeros(1,length(CC));
            for ci = 1:length(CC)
                cvx_begin
                    variables teqTop(15) teq(N) lam(lamdim)
                    minimize( CC(ci)*norm(lam,1)+50*norm(teq,2)+100*norm(teqTop,2))
                    subject to
                    Aeq*lam==beq+teq
                    AeqTop*lam==beqTop+teqTop
                    Aineq*lam <= bineq
                cvx_end
                
                LAMS(:,ci)=lam;
                costs(ci) = norm(teq,2);
            end
            [~,bind] = min(costs);
            lam = LAMS(:,bind);
            lam(1) = lam(1)+log(factconst);
            lamsol = lam;
            
            mxentpoly=zeros(1,dim+1);
            for i=1:lamdim
                mxentpoly=add_sub_polyND(mxentpoly, scalar_multiply_polyND(lamsol(i),Pf{i}),'add');
            end
            mxentpoly=simplify_polyND(mxentpoly);
            
            obj.mxentpoly = mxentpoly;

        end
        
        function PlotExpPolyFits(obj,states,LB,UB)
                fest=evaluate_polyND(obj.mxentpoly,obj.X);
                ftrue = log(obj.p);
                
                plot3(obj.X(:,states(1)),obj.X(:,states(2)),ftrue,'ro' )
                hold on
                plot3(obj.X(:,states(1)),obj.X(:,states(2)),fest,'bs' )
                for i =1:3
                    Xmc = mvurnd(LB,UB,5000);
                    fmc=evaluate_polyND(obj.mxentpoly,Xmc);
                    indmcbad = fmc>max(ftrue);
                    fmc(indmcbad);
                    plot3(Xmc(indmcbad,states(1)),Xmc(indmcbad,states(2)),fmc(indmcbad),'g+' )
                end
        end
    end
end

% Xtrue and 
