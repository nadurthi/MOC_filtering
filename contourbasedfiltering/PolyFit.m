
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
        
        function mxentpoly = fitExpPoly_A_Atop_Aineq(obj,Pf,Xineq,Xtest)
            [N,dim] = size(obj.X);
            lamdim = length(Pf);
            Aeq=zeros(N,lamdim);
            
%             for r=1:1:N
%                 Aeq(r,:) = evaluate_MatrixOfPolys(Pf,obj.X(r,:));
%             end
            tic
            for ib=1:length(Pf)
                Aeq(:,ib)=evaluate_PolyforLargetSetX(Pf{ib},obj.X);
            end
            toc
            
            
            factconst = max(obj.p)/10;
            pnfit = obj.p/factconst;
            beq = log(pnfit);
            
            Ntop = 15;
            AeqTop = Aeq(1:Ntop,:);
            beqTop = beq(1:Ntop);
            
            Atest=zeros(size(Xtest,1),lamdim);
            for ib=1:length(Pf)
                Atest(:,ib)=evaluate_PolyforLargetSetX(Pf{ib},Xtest);
            end
            btest = max(beq)*ones(size(Atest,1),1);
            
            
            % lamdim=length(lam);
%             keyboard
            Nineq = size(Xineq,1);
            Aineq=zeros(size(Xineq,1),lamdim);
            tic
            for ib=1:length(Pf)
                Aineq(:,ib)=evaluate_PolyforLargetSetX(Pf{ib},Xineq);
            end
            toc
%             tic
%             for r=1:1:size(Xineq,1)
%                 Aineq(r,:) = evaluate_MatrixOfPolys(Pf,Xineq(r,:));
%             end
%             toc
            bineq = 1*(min(beq)-0 )*ones(size(Aineq,1),1);
            
            
            
            % %working good
            %     minimize( 10*norm(lam2,1)+50*norm(t,2)+150*norm(t2,2))
            CC=linspace(0,100,21);
            LAMS=zeros(lamdim,length(CC));
            costs = zeros(1,length(CC));
            for ci = 1:length(CC)
                cvx_begin
                    variables teq(N) lam(lamdim)
                    minimize( CC(ci)*norm(lam,1)+100*norm(teq,2) )
                    subject to
                    Aeq*lam==beq+teq;
%                     AeqTop*lam==beqTop+teqTop;
                    Aineq*lam <=bineq;
                cvx_end
                
                if all(Atest*lam <= btest)
                    disp('All probs are in the bounds')
                    break
                end

                LAMS(:,ci)=lam;
                costs(ci) = norm(teq,2);
            end
%             [~,bind] = min(costs);
%             lam = LAMS(:,bind);
            lam(1) = lam(1)+log(factconst);
            lamsol = lam;
            
            mxentpoly=zeros(1,dim+1);
            for i=1:lamdim
                mxentpoly=add_sub_polyND(mxentpoly, scalar_multiply_polyND(lamsol(i),Pf{i}),'add');
            end
            mxentpoly=simplify_polyND(mxentpoly);
            
            obj.mxentpoly = mxentpoly;

        end
        function mxentpoly = fitExpPoly_A_Atop_AdaptiveC(obj,Pf,Xtest)
            [N,dim] = size(obj.X);
            lamdim = length(Pf);
            Aeq=zeros(N,lamdim);

            for ib=1:length(Pf)
                Aeq(:,ib)=evaluate_PolyforLargetSetX(Pf{ib},obj.X);
            end

            factconst = max(obj.p)/10;
            pnfit = obj.p/factconst;
            beq = log(pnfit);
            
            Ntop = 15;
            AeqTop = Aeq(1:Ntop,:);
            beqTop = beq(1:Ntop);

            Atest=zeros(size(Xtest,1),lamdim);
            for ib=1:length(Pf)
                Atest(:,ib)=evaluate_PolyforLargetSetX(Pf{ib},Xtest);
            end
            btest = max(beq)*ones(size(Atest,1),1);

            CC=linspace(0,100,21);
            LAMS=[];
            costs = [];
            for ci = 1:length(CC)
                cvx_begin
                    variables teq(N) lam(lamdim) teqTop(15) 
                    minimize( CC(ci)*norm(lam,1)+100*norm(teq,2)) 
                    subject to
                    Aeq*lam==beq+teq;
                    AeqTop*lam==beqTop+teqTop;
                    AeqTop(1,:)*lam==beqTop(1);
                cvx_end
                
                LAMS = [LAMS,lam(:)];
                costs = [costs,norm(teq,2)];
                if all(Atest*lam <= btest)
                    disp('All probs are in the bounds')
                    break
                end
                
%                 keyboard
            end
            lamsol = lam;
            
            mxentpoly=zeros(1,dim+1);
            for i=1:lamdim
                mxentpoly=add_sub_polyND(mxentpoly, scalar_multiply_polyND(lamsol(i),Pf{i}),'add');
            end
            mxentpoly=simplify_polyND(mxentpoly);
            
%             mxentpoly2=zeros(length(Pf),dim+1);
%             mxentpoly2(:,1)=lamsol;
%             for i=1:lamdim
%                 mxentpoly2(i,2:end) = Pf{i}(2:end);
%             end
            c0 = get_coeff_NDpoly(mxentpoly,zeros(1,dim));
            mxentpoly = update_or_insert_coeff_NDpoly(mxentpoly,zeros(1,dim),c0+log(factconst));
            
            
%             keyboard
            
%             fest=evaluate_polyND(mxentpoly,obj.X);
% %             fest2=evaluate_polyND(mxentpoly2,obj.X);
%             figure
%             plot3(obj.X(:,1),obj.X(:,2),beq,'ro')
%             hold on
%             plot3(obj.X(:,1),obj.X(:,2),fest,'bs')
% %             plot3(obj.X(:,1),obj.X(:,2),fest2,'g+')
            
            
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
        function PlotExpPolyFits_points(obj,states,LB,UB)
                fest=evaluate_polyND(obj.mxentpoly,obj.X);
                ftrue = log(obj.p);
                
                plot3(obj.X(:,states(1)),obj.X(:,states(2)),ftrue,'ro' )
                hold on
                plot3(obj.X(:,states(1)),obj.X(:,states(2)),fest,'bs' )
                for i =1:3
                    Xmc = mvurnd(LB,UB,5000);
                    fmc=evaluate_polyND(obj.mxentpoly,Xmc);
%                     indmcbad = fmc>max(ftrue);
%                     fmc(indmcbad);
                    plot3(Xmc(:,states(1)),Xmc(:,states(2)),fmc,'g+' )
                end
        end
        
    end
end

% Xtrue and 
