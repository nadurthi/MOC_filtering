classdef DataSet < handle
    %maintain the points and their probabilities
    % provide methods to trasnform them
    
    properties
        X;
        p;
        Xorig;
        porig;
        affinetransforms={};
        cntrans=0;
        originalState='TrueState';
        sortind;
    end
    methods
        % constructor. save a copy constants
        function obj= DataSet(X,p,originalState)
            % 'TrueState' or 'TransformedState'
            obj.X=X;
            obj.p=p;
            obj.originalState=originalState;
            
            obj.Xorig=X;
            obj.porig=p;
        end
        function ds = MakeCopy(obj)
            ds = DataSet(obj.X,obj.p,obj.originalState);
            ds.affinetransforms=obj.affinetransforms;
            ds.cntrans=0;
            ds.originalState='TrueState';
        end
        function [mm,PP] = GetMeanCov(obj)
           [mm,PP] = MeanCov(obj.X,obj.p/sum(obj.p)); 
        end
        function obj = AddAffineTrasform(obj,A,m)
            m=m(:);
            obj.cntrans=obj.cntrans+1;
            obj.affinetransforms{obj.cntrans}={A,m};
            
            for i=1:size(obj.X,1)
                obj.X(i,:) = A*obj.X(i,:)'+m(:);
            end
            
            obj.p =  obj.p/det(A);
        end
        function obj = AddMeanCov_to_OI_Trasform(obj,mx,Px)
            % remember the P is cov of X nad sqrtm(inv(P)) is used
%             A = inv(sqrtm(Px));
            [u,v]=eig(Px);
            vv=1./diag(sqrt(v));
            A=u*diag(vv)*u';
            
            m = -A*mx;
            obj = obj.AddAffineTrasform(A,m);
        end
        
        function obj = AddHyperCubeTrasform(obj,LB,UB)
            maxX = max(obj.X,[],1);
            minX = min(obj.X,[],1);
            maxX=maxX(:);
            minX=minX(:);
            
            LB=LB(:);
            UB=UB(:);
            A = diag((UB-LB)./(maxX-minX));
            m = -(minX.*(UB-LB))./(maxX-minX)+LB;
            
            obj = obj.AddAffineTrasform(A,m);
            
        end
        function [A,m]=GetAffineTransform_Original2Final(obj)
            % y =A4 A3 A2 A1 x + A4 A3 A2 m1 + A4 A3 m2 + A4 m3 + m4
            A=1;
            for i=1:obj.cntrans
                A=obj.affinetransforms{i}{1}*A;
            end
            m=0;
            for i=obj.cntrans:-1:1
                if i==obj.cntrans
                    c=1;
                else
                    c=c*obj.affinetransforms{i+1}{1};
                end
                m=m+c*obj.affinetransforms{i}{2};
            end
        end
        function [A,m]=GetAffineTransform_Final2Original(obj)
            [Aof,mof]=obj.GetAffineTransform_Original2Final();
            A=inv(Aof);
            m=-A*mof;
        end
        function transForms = GetTrasnformers(obj)
            [Atranf2norm,mtransf2norm]=obj.GetAffineTransform_Original2Final();
            [Atranf2true,mtransf2true]=obj.GetAffineTransform_Final2Original();
            transForms.A_true2norm = Atranf2norm;
            transForms.m_true2norm = mtransf2norm;
            transForms.A_norm2true = Atranf2true;
            transForms.m_norm2true = mtransf2true;
            
            transForms.trueX2normX = @(x)affineTransform(x,Atranf2norm,mtransf2norm);
            transForms.normX2trueX = @(xn)affineTransform(xn,Atranf2true,mtransf2true);
            transForms.normprob2trueprob = @(p)p/det(Atranf2true);
            transForms.trueprob2normprob = @(p)p/det(Atranf2norm);
            
            transForms.normprob2trueprob_constant = 1/det(Atranf2true);
            transForms.trueprob2normprob_constant = 1/det(Atranf2norm);
            
        end
        function obj = Reset2Original(obj)
            [A,m]=obj.GetAffineTransform_Final2Original();
            for i=1:size(obj.X,1)
                obj.X(i,:) = A*obj.X(i,:)'+m(:);
            end
            
            obj.p =  obj.p/det(A);
            
            obj.affinetransforms={};
            obj.cntrans=0;
        end
        function [Y,py] = GetOriginal(obj)
            % if there are no transforms just return X
            if isempty(obj.affinetransforms)
                Y=obj.X;
                py=obj.p;
                
                return;
            end
            
            % here assuming that original X was trasnformed
            [A,m]=obj.GetAffineTransform_Final2Original();
            Y=zeros(size(obj.X));
            
            for i=1:size(obj.X,1)
                Y(i,:) = A*obj.X(i,:)'+m(:);
            end
            
            py =  obj.p/det(A);
            
        end
        function [Yt,pt] = ApplyAffineTransform_Original2Final(obj,Y,p)
            [r,c] = size(Y);
            if r==1 || c==1
                Y=Y(:)';
            end
            p=p(:);
            
            Yt=zeros(size(Y));
            
            [A,m]=obj.GetAffineTransform_Original2Final();
            
            for i=1:size(Y,1)
                Yt(i,:) = A*Y(i,:)'+m(:);
            end
            
            pt =  p/det(A);
        end
        function [Yt,pt] = ApplyAffineTransform_Final2Original(obj,Y,p)
            [r,c] = size(Y);
            if r==1 || c==1
                Y=Y(:)';
            end
            p=p(:);
            
            Yt=zeros(size(Y));
            
            [A,m]=obj.GetAffineTransform_Final2Original();
            
            for i=1:size(Y,1)
                Yt(i,:) = A*Y(i,:)'+m(:);
            end
            
            pt =  p/det(A);
        end
        function pt = ApplyAffineTransformProb_Original2Final(obj,p)
            
            p=p(:);
            
            [A,m]=obj.GetAffineTransform_Original2Final();
            
            pt =  p/det(A);
        end
        function pt = ApplyAffineTransformProb_Final2Original(obj,p)
            
            p=p(:);
            
            [A,m]=obj.GetAffineTransform_Final2Original();
            
            pt =  p/det(A);
        end
        
        % Some helpful functions
        function obj = SortByProb(obj,dir)
            [~,ind]=sort(obj.p(:),1,dir);
            obj.sortind=ind;
            obj.p=obj.p(ind);
            obj.X=obj.X(ind,:);
            
        end
        function logp = getlogProb(obj)
            logp = log(obj.p);
        end
        function obj = PlotPoints2D(obj,states,c)
            plot(obj.X(:,states(1)),obj.X(:,states(2)),c)
        end
        function obj = PlotPoints3D(obj,states,c)
            plot3(obj.X(:,states(1)),obj.X(:,states(2)),obj.X(:,states(3)),c)
        end
        function obj = PlotPointsProbs3D(obj,states,c)
            plot3(obj.X(:,states(1)),obj.X(:,states(2)),obj.p,c)
        end
        
        
    end
end



