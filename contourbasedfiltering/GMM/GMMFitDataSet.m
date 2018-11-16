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
        
        function GMM = FitGMM_kmeans_optimwt(obj,NgcompMax)
            % optimize the wts of GMM(computed from FitGMM_kmean_equalwt) to the probabilituy
            [Np,dim]=size(obj.X);
            GMM=obj.FitGMM_kmean_equalwt(NgcompMax);
            Ngcomp = GMM.Ngcomp;
            Ag = zeros(Np,Ngcomp);
            for i=1:Ngcomp
                Ag(:,i) = mvnpdf(obj.X,GMM.mx{i}',GMM.Px{i});
            end
            bg = obj.p;
            
            cvx_begin
            variables wg(Ngcomp) t(Np)
            minimize( norm(Ag*wg-bg,1) )
            subject to
            wg>=0;
            sum(wg)==1;
            cvx_end
            wg = wg+0.0001;
            wg = wg/sum(wg);
            GMM.w = wg(:);
            obj.GMM = GMM;
        end
        %%
        function GMM = FitGMM_BoostingGaussian(obj,Ngcomp,ErrWt)
            GMM = BoostFit(obj.X,obj.p,Ngcomp);
            obj.GMM = GMM;
            
        end
        %% Hulling operations
        function GMMhull = SetGMM_Hull(obj,NgcompMax)
            obj.GMMhull = GetGMM_Hull_impl(obj.X,NgcompMax,'kmeans');
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
            lb1 = mm(states(1)) - 3*PPsqrt(states(1),states(1));
            ub1 = mm(states(1)) + 3*PPsqrt(states(1),states(1));
            
            lb2 = mm(states(2)) - 3*PPsqrt(states(2),states(2));
            ub2 = mm(states(2)) + 3*PPsqrt(states(2),states(2));
            
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
            plot3(obj.X(:,states(1)),obj.X(:,states(2)),obj.p,c,'MarkerSize',7)
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




