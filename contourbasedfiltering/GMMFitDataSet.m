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
%             flg=1;
%             prevIDX=0;
%             for Nclust = 1:NgcompMax
%                 [IDX,C] = kmeans(obj.X, Nclust);
%                 for i=1:size(C,1)
%                     [m,pR]=MeanCov(obj.X(IDX==i,:),obj.p(IDX==i)/sum(obj.p(IDX==i)));
%                     egs = eig(pR);
%                     if any(egs<0) || cond(pR)>5000
%                         flg=0;
%                         break
%                     end
%                     if flg==0
%                         break;
%                     end
%                 end
%                 if flg==0
%                     break;
%                 end
%                 prevIDX = IDX;
%                 
%             end
            IDX = GenerateClusterIndexes(obj.X,NgcompMax,'gmm');
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
            minimize( norm(Ag*wg-bg) )
            subject to
            wg>=0
            sum(wg)==1
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
        function GMMhull = SetGMM_Hull(obj)
            NgcompMax = 3;
            %             obj.GMMhull = GetGMM_kmeans_Hull_impl(obj.X);
%             obj.GMMhull = GetGMM_AggClust_Hull_impl(obj.X);
            obj.GMMhull = GetGMM_Hull_impl(obj.X,NgcompMax,'gmm');
            GMMhull = obj.GMMhull;
        end
        
        function ind = IsInsideHull(obj,Xtest,factor)
            ind = zeros(size(Xtest,1),1);
            for nc = 1:obj.GMMhull.Ngcomp
                mm = obj.GMMhull.mx{nc};
                PP = obj.GMMhull.Px{nc};
                ind = ind | CheckifInsideEllipsoid(Xtest,mm,factor*PP);
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
                        GMMprobs(i,j) = GMMprobs(i,j) + obj.GMM.w(nk)*mvnpdf([xx(i,j),yy(i,j)],obj.GMM.mx{nk}',obj.GMM.Px{nk});
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
        
        
        function plotGMMpointsHUll(obj,states,factor,c)
            plot(obj.X(:,states(1)),obj.X(:,states(2)),c)
            hold on
            for i=1:obj.GMMhull.Ngcomp
                plot_nsigellip(obj.GMMhull.mx{i}(states),factor*obj.GMMhull.Px{i}(states,states),1,'k',1)
            end
        end
        
        
        
        
    end
end

function GMM = BoostFit(X,p,Ngcomp)

[N,dim] = size(X);

end

function idx = GenerateClusterIndexes(X,NgcompMax,method)
[N,dim] = size(X);


previdx=[];
for Ngcomp=3:NgcompMax
    
    if strcmp(method,"kmeans")
        idx = kmeans(X,Ngcomp);
    elseif strcmp(method,"AggClust")
        Z = linkage(X,'complete');
        idx = cluster(Z,'Maxclust',Ngcomp);
    elseif strcmp(method,"gmm")
        GMModel = fitgmdist(X,Ngcomp);
        idx = cluster(GMModel,X) ;
    end
    
    flg=1;
    for i=1:Ngcomp
        xx=X(idx==i,:) ;
        Nc = size(xx,1);
        if Nc <= dim
            disp('BREAK: Nc <= dim')
            flg=0;
            break
        end
        ww = ones(Nc,1)/Nc;
        [mcp,pcp]=MeanCov(xx,ww/sum(ww));
        eigsP = eig(pcp);
        if min(eigsP) <=0 || ~isreal(eigsP) || any(isnan(eigsP)==1) || any(isfinite(eigsP)==0)
            disp('BREAK: min(eigsP) <=0 || isreal(eigsP) || any(isnan(eigsP)==1) || any(isfinite(eigsP)==0)')
            flg=0;
            break
        end
        if cond(pcp)>5000
            disp('***********\n cond(pcp)>5000 \n *******************')
            disp(Ngcomp)
            flg=0;
            break
        end
    end
    if isempty(previdx)
        previdx = idx;
    end
    if flg==1
        previdx = idx;
    else
        break
    end
    
end

idx = previdx;

end

function GMMhull = GetGMM_Hull_impl(X,NgcompMax,method)
[N,dim] = size(X);

idx = GenerateClusterIndexes(X,NgcompMax,method);
Ngcomp = max(idx);
MF =cell(Ngcomp,1);
PF =cell(Ngcomp,1);

for i=1:Ngcomp
    xx=X(idx==i,:) ;
    Nc = size(xx,1);
    ww = ones(Nc,1)/Nc;
    [mcp,pcp]=MeanCov(xx,ww/sum(ww));
    xx = xx';
    [n,m] = size(xx);
    cvx_begin
        variable A(n,n) symmetric
        variable b(n)
        maximize( det_rootn( A ) )
        subject to
            norms( A * xx + b * ones( 1, m ), 2 ) <= 1;
    cvx_end

%     for c = 1:1:1000
%         ind = CheckifInsideEllipsoid(xx,mcp,c*pcp);
%         if prod(ind)==1
%             break
%         end
%     end
% keyboard
 
    MF{i} = -inv(A)*b;
    PF{i} = inv(A);
    
end
GMMhull.w = ones(Ngcomp,1)/Ngcomp;
GMMhull.mx = MF;
GMMhull.Px = PF;
GMMhull.Ngcomp = Ngcomp;

end

% 
% function GMMhull = GetGMM_kmeans_Hull_impl(X)
% [N,dim] = size(X);
% MF ={};
% PF ={};
% 
% for Ngcomp=1:10
%     
%     idx = kmeans(X,Ngcomp);
%     Mcomps =cell(Ngcomp,1);
%     Pcomps =cell(Ngcomp,1);
%     Ag=zeros(N,Ngcomp);
%     flg=1;
%     for i=1:Ngcomp
%         xx=X(idx==i,:) ;
%         Nc = size(xx,1);
%         if Nc <= dim
%             disp('BREAK: Nc <= dim')
%             flg=0;
%             break
%         end
%         ww = ones(Nc,1)/Nc;
%         [mcp,pcp]=MeanCov(xx,ww/sum(ww));
%         eigsP = eig(pcp);
%         if min(eigsP) <=0 || ~isreal(eigsP) || any(isnan(eigsP)==1) || any(isfinite(eigsP)==0)
%             disp('BREAK: min(eigsP) <=0 || isreal(eigsP) || any(isnan(eigsP)==1) || any(isfinite(eigsP)==0)')
%             flg=0;
%             break
%         end
%         if cond(pcp)>1000
%             disp('cond(pcp)>1000')
%             flg=0;
%             break
%         end
%         for c = 1:1:100
%             ind = CheckifInsideEllipsoid(xx,mcp,c*pcp);
%             if prod(ind)==1
%                 break
%             end
%         end
%         Mcomps{i} = mcp(:);
%         Pcomps{i} = c*pcp;
%     end
%     if flg==1
%         MF = Mcomps;
%         PF = Pcomps;
%     else
%         break
%     end
% end
% 
% GMMhull.w = ones(length(MF),1)/length(MF);
% GMMhull.mx = MF;
% GMMhull.Px = PF;
% GMMhull.Ngcomp = length(MF);
% end

function ind = CheckifInsideEllipsoid(X,m,P)
m=m(:);
invP  = inv(P);
const = 1/sqrt(det(2*pi*P));
[N,dim] = size(X);

xt = sqrtm(P)*[1;zeros(dim-1,1)]+m;
pt = const*exp(-0.5*(xt-m)'*invP*(xt-m));

p = zeros(N,1);
for i=1:N
    x = X(i,:)';
    p(i) = const*exp(-0.5*(x-m)'*invP*(x-m));
end

ind = zeros(N,1);
ind(p>pt)=1;

end





