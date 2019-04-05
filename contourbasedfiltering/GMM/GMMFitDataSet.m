classdef GMMFitDataSet < handle
    % fit GMMs and plot the visualizations
    
    properties
        X;
        p;
        GMM;
        GMMhull;
    end
    methods
        % constructor. save a copy constants
        function obj= GMMFitDataSet(X,p)
            % 'TrueState' or 'TransformedState'
            obj.X=X;
            obj.p=p;
        end
        
        function GMM = FitGMM_1comp(obj)
            % DEPRECATED: Specify use case before using as it si notc clear
            [N,dim]=size(obj.X);
            Ntop = max(floor(N/10),2*dim+5);
            [mm,PP] = MeanCov(obj.X(1:Ntop,:),obj.p(1:Ntop)/sum(obj.p(1:Ntop)));
            c=fmincon(@(c)max( (obj.p(1:Ntop)-mvnpdf(obj.X(1:Ntop,:),mm(:)',c*PP))./obj.p(1:Ntop)),1,[],[],[],[],0.0001,1 );
            PP = c*PP;
            for w=10:-0.01:0.01
                if all(obj.p>w*mvnpdf(obj.X,mm(:)',PP))
                    break
                end
            end
            [c,w,all(obj.p>w*mvnpdf(obj.X,mm(:)',PP))]
            obj.GMM.w = w;
            obj.GMM.mx = {mm};
            obj.GMM.Px = {PP};
            obj.GMM.Ngcomp = 1;
            
            GMM = obj.GMM;
        end
        %%
        function GMM = FitGMM_kmean_equalwt(obj,NgcompMax)
            % Do kmeans and just compute the mean and covariance of  those
            % cluster
            IDX = GenerateClusterIndexes(obj.X,NgcompMax,'kmeans');
            %             IDX = prevIDX;
            Nclust = max(IDX);
            obj.GMM.w = ones(Nclust,1)/Nclust;
            obj.GMM.mx = cell(Nclust,1);
            obj.GMM.Px = cell(Nclust,1);
            obj.GMM.Ngcomp = Nclust;
            
            for i=1:Nclust
                [m,pR]=MeanCov(obj.X(IDX==i,:),obj.p(IDX==i)/sum(obj.p(IDX==i)));
                obj.GMM.mx{i} = m;
                obj.GMM.Px{i} = pR;
            end
            
            GMM = obj.GMM;
        end
        
        function GMM = FitGMM_dummy(obj,X,NgcompMax,mineigvalset)
            % optimize the wts of GMM(computed from FitGMM_kmean_equalwt) to the probabilituy
            idx = GenerateClusterIndexes_merging(X,NgcompMax,'None',mineigvalset);
            Ngcomp = max(idx);
            
            MF =cell(Ngcomp,1);
            PF =cell(Ngcomp,1);
            AF=cell(Ngcomp,1);
            BF=cell(Ngcomp,1);

            for i=1:Ngcomp
                xx=X(idx==i,:) ;
                Nc = size(xx,1);
                ww = ones(Nc,1)/Nc;
                [mcp,pcp]=MeanCov(xx,ww/sum(ww));
                MF{i} = mcp;
                PF{i} = pcp;
                [u,s]=eig( PF{i} );
                sd = diag(s);
                sd(sd<0.01)=0.01;
                s=diag(sd);
                PF{i} = u*s*u';
                A =sqrtm(inv(PF{i}));
                b=-A*MF{i};
                AF{i} = A;
                BF{i} = b;
            end
            GMM.w = ones(Ngcomp,1)/Ngcomp;
            GMM.mx = MF;
            GMM.Px = PF;
            GMM.Ngcomp = Ngcomp;
            GMM.A = AF;
            GMM.b = BF;
            GMM.idx=idx;

            obj.GMM=GMM;
                
        end
        %%
        function GMM = FitGMM_BoostingGaussian(obj,Ngcomp,ErrWt)
            GMM = BoostFit(obj.X,obj.p,Ngcomp);
            obj.GMM = GMM;
            
        end
        %% Hulling operations
        function GMMhull = SetGMM_Hull(obj,NgcompMax,mineigvalsetsqr)
            switch nargin
                case 3
                    mineigvalsqr=mineigvalsetsqr;
                otherwise
                    mineigvalsqr=0.1;
            end

            obj.GMMhull = GetGMM_Hull_impl(obj.X,NgcompMax,'kmeans',mineigvalsqr);
            GMMhull = obj.GMMhull;
        end
        function obj = resetGMM_Hull_wts(obj)
            obj.GMMhull.w=obj.GMMhull.w/sum(obj.GMMhull.w);
        end
        function obj = resetGMM_Hull_full(obj)
            obj.GMMhull.w=[];
            obj.GMMhull.mx={};
            obj.GMMhull.Px={};
            obj.GMMhull.Ngcomp=0;
            obj.GMMhull.A={};
            obj.GMMhull.b={};
        end
        function GMMhull = optimGMMhullwts(obj,pp)
            switch nargin
                
                case 2
                    bg=pp;
                case 1
                    bg = obj.p;
                otherwise
                    error('aruments for GMMhullwt optimization not valid');
            end

            [Np,dim]=size(obj.X);
            Ngcomp = obj.GMMhull.Ngcomp;
            Ag = zeros(Np,Ngcomp);
            for i=1:Ngcomp
                Ag(:,i) = mvnpdf(obj.X,obj.GMMhull.mx{i}',obj.GMMhull.Px{i});
            end
%             bg = obj.p;
            
            cvx_begin
                variables wg(Ngcomp) t(Np)
                minimize( norm((Ag*wg-bg)./bg,4) )
                subject to
                wg>=0;
                sum(wg)==1;
            cvx_end
            wg = wg+0.0001;
            wg = wg/sum(wg);
            obj.GMMhull.w = wg(:);

            
            GMMhull = obj.GMMhull;
        end
        function GMMhull = optimGMMhullwts_relative2probs(obj,XX,pp)
            w=obj.GMMhull.w*0;
            for i=1:obj.GMMhull.Ngcomp
                A=obj.GMMhull.A{i};
                b=obj.GMMhull.b{i};
                ind = CheckifInsideEllipsoid_Abmethod(XX,A,b,1);
                w(i)=max(pp(ind==1));
            end
            w = w+0.0001;
            w = w/sum(w);
            obj.GMMhull.w = w(:);            
            GMMhull = obj.GMMhull;
            
        end
        function GMM = optimGMMwts_relative2probs_and_setmXPx(obj,XX,pp)
            w=obj.GMM.w*0;
            for i=1:obj.GMM.Ngcomp
                A=obj.GMM.A{i};
                b=obj.GMM.b{i};
                ind = CheckifInsideEllipsoid_Abmethod(XX,A,b,1);
                xx=XX(ind==1,:);
                pq=pp(ind==1);
                [mX,PX]=MeanCov(xx,pq/sum(pq));
                obj.GMM.mx{i}=mX;
                obj.GMM.Px{i}=2*PX;
                w(i)=max(pq);
            end
            w = w+0.0001;
            w = w/sum(w);
            obj.GMM.w = w(:);            
            GMM = obj.GMM;
            
        end
        function GMM = optimGMMwts_relative2probs(obj,XX,pp)
            w=obj.GMM.w*0;
            for i=1:obj.GMM.Ngcomp
                A=obj.GMM.A{i};
                b=obj.GMM.b{i};
                ind = CheckifInsideEllipsoid_Abmethod(XX,A,b,1);
                pq=pp(ind==1);
                w(i)=mean(pq);
            end
            w = w+0.0001;
            w = w/sum(w);
            obj.GMM.w = w(:);            
            GMM = obj.GMM;
            
        end
        
        function GMMhull = optimGMMhullwts_reoptimize(obj,XX,pp)


            [Np,dim]=size(XX);
            Ngcomp = obj.GMMhull.Ngcomp;
            Ag = zeros(Np,Ngcomp);
            for i=1:Ngcomp
                Ag(:,i) = mvnpdf(XX,obj.GMMhull.mx{i}',obj.GMMhull.Px{i});
            end
            bg = pp;
%             wg=fmincon(@(wg)norm((Ag*wg-bg)./(bg+1),4),ones(Ngcomp,1)/Ngcomp,[],[],ones(1,Ngcomp),1,zeros(Ngcomp,1),ones(Ngcomp,1) );
            cvx_begin
                variables wg(Ngcomp) t(Np)
                minimize( norm((Ag*wg-bg)./(bg+1),2) )
                subject to
                wg>=0;
%                 sum(wg)==1;
            cvx_end
            wg = wg+0.0001;
            wg = wg/sum(wg);
            obj.GMMhull.w = wg(:);

            
            GMMhull = obj.GMMhull;
        end
        
        function ind = IsInsideHull(obj,Xtest,factor)
            ind = zeros(size(Xtest,1),1);
            %             for nc = 1:obj.GMMhull.Ngcomp
            %                 mm = obj.GMMhull.mx{nc};
            %                 PP = obj.GMMhull.Px{nc};
            %                 ind = ind | CheckifInsideEllipsoid(Xtest,mm,factor*PP);
            %             end
            for nc = 1:obj.GMMhull.Ngcomp
                A = obj.GMMhull.A{nc};
                b = obj.GMMhull.b{nc};
                ind = ind | CheckifInsideEllipsoid_Abmethod(Xtest,A,b,factor);
            end
            
            
        end
        %% sampling the GMM
        function X = gen_quad_GMM_GH_2D(obj,Nqmax)
%             switch(method)
%                 case 'UT'
%                     qd_pts=UT_sigmapoints(m,P,2);
%                 case 'GH'
%                     qd_pts=@(m,P)GH_pts(m,P,para);
%             end
            X=[];
            for i=1:obj.GMM.Ngcomp
               if obj.GMM.w(i)>0.7
                   Npara=9;
                   c=0.2;
               elseif obj.GMM.w(i)>0.5 && obj.GMM.w(i)<=0.7  
                   Npara=9;
                   c=0.2;
               elseif obj.GMM.w(i)>0.3 && obj.GMM.w(i)<=0.5
                   Npara=9;
                   c=0.01;
               elseif obj.GMM.w(i)>0.1 && obj.GMM.w(i)<=0.3
                   Npara=4;
                   c=0.4;
               elseif obj.GMM.w(i)>1e-2 && obj.GMM.w(i)<=0.1
                   Npara=3;
                   c=0.4;
               elseif obj.GMM.w(i)>1e-4 && obj.GMM.w(i)<=1e-2
                   Npara=2;
                   c=0.4;
               elseif obj.GMM.w(i)<1e-4
                   Npara=0;
                   c=0.4;
               end
               if Npara==0
                   x=[];
               else
                    x = GH_pts(obj.GMM.mx{i},c^2*obj.GMM.Px{i},Npara);
               end
               X=vertcat(X,x);
            end
      
        end
        function X = gen_quad_GMM_GH_4D(obj,Nqmax)
%             switch(method)
%                 case 'UT'
%                     qd_pts=UT_sigmapoints(m,P,2);
%                 case 'GH'
%                     qd_pts=@(m,P)GH_pts(m,P,para);
%             end
            X=[];
            for i=1:obj.GMM.Ngcomp
               if obj.GMM.w(i)>0.5
                   Npara=9;
                   c=0.2;
               elseif obj.GMM.w(i)>0.3 && obj.GMM.w(i)<=0.5  
                   Npara=9;
                   c=0.2;
               elseif obj.GMM.w(i)>0.1 && obj.GMM.w(i)<=0.3
                   Npara=8;
                   c=0.2;
               elseif obj.GMM.w(i)>1e-2 && obj.GMM.w(i)<=0.1
                   Npara=5;
                   c=0.2;
               elseif obj.GMM.w(i)>1e-5 && obj.GMM.w(i)<=1e-2
                   Npara=3;
                   c=0.2;
               elseif obj.GMM.w(i)<1e-5
                   Npara=0;
                   c=0.2;
               end
               if Npara==0
                   x=[];
               else
                    x = GH_pts(obj.GMM.mx{i},c^2*obj.GMM.Px{i},Npara);
               end
               X=vertcat(X,x);
            end
      
        end
        %%
        
        function plotGMMpoints(obj,states,c)
            plot(obj.X(:,states(1)),obj.X(:,states(2)),c)
            hold on
            for i=1:obj.GMM.Ngcomp
                plot_nsigellip(obj.GMM.mx{i}(states),obj.GMM.Px{i}(states,states),1,'r',1)
            end
        end
        function plotGMMSurf(obj,states,c)
            [mm,PP] = MeanCov(obj.X,obj.p/sum(obj.p));
            PPsqrt = sqrtm(PP);
            lb1 = mm(states(1)) - 5*PPsqrt(states(1),states(1));
            ub1 = mm(states(1)) + 5*PPsqrt(states(1),states(1));
            
            lb2 = mm(states(2)) - 5*PPsqrt(states(2),states(2));
            ub2 = mm(states(2)) + 5*PPsqrt(states(2),states(2));
            
            [xx,yy] = meshgrid(linspace(lb1,ub1,70 ),linspace(lb2,ub2,70 ) );
            GMMprobs = zeros(size(xx));
            for i=1:size(xx,1)
                for j=1:size(xx,2)
                    for nk=1:obj.GMM.Ngcomp
                        GMMprobs(i,j) = GMMprobs(i,j) + obj.GMM.w(nk)*mvnpdf([xx(i,j),yy(i,j)],obj.GMM.mx{nk}(states)',obj.GMM.Px{nk}(states,states));
                    end
                end
            end
            mesh(xx,yy,GMMprobs);
            alpha 0.5
            hold on;
            %             plot(obj.X(:,states(1)),obj.X(:,states(2)),c)
%             plot3(obj.X(:,states(1)),obj.X(:,states(2)),obj.p,c,'MarkerSize',7)
            for i=1:obj.GMM.Ngcomp
                plot_nsigellip(obj.GMM.mx{i}(states),obj.GMM.Px{i}(states,states),1,'r',1)
            end
        end
        
        
        function plotGMMpointsHUll(obj,states,Xineq,factor,c)
            plot(obj.X(:,states(1)),obj.X(:,states(2)),c)
            hold on
            for i=1:obj.GMMhull.Ngcomp
                plot_nsigellip(obj.GMMhull.mx{i}(states),factor*obj.GMMhull.Px{i}(states,states),1,'k',1)
            end
             if isempty(Xineq)==0
                    plot(Xineq(:,states(1)),Xineq(:,states(2)),'b+')
             end
             hold off
        end
        
        function plotGMMhullPoints_debug(obj,states,Xineq,factor,c)
            
            
            for i=1:obj.GMMhull.Ngcomp
                y = obj.GMMhull.idx==i;
                plot(obj.X(y,states(1)),obj.X(y,states(2)),c)
                hold on
                if isempty(Xineq)==0
                    plot(Xineq(:,states(1)),Xineq(:,states(2)),'b+')
                end
                plot_nsigellip(obj.GMMhull.mx{i}(states),factor*obj.GMMhull.Px{i}(states,states),1,'r',1)
                hold off
                
                axis([-2,2,-2,2])
                pause(0.5)
            end
        end
%         
        
    end
end




