classdef GMMFitDataSet < handle
    % fit GMMs and plot the visualizations
    
    properties
        X;
        p;
    end
    methods
        % constructor. save a copy constants
        function obj= GMMFitDataSet(X,p)
            % 'TrueState' or 'TransformedState'
            obj.X=X;
            obj.p=p;
        end
        
        function GMM = FitGMM_kmeans_NoWeightConstraint(obj,Ngcomp,ErrWt)
            [N,dim] = size(obj.X);
            idx = kmeans(obj.X,Ngcomp);
            Mcomps =cell(Ngcomp,1);
            Pcomps =cell(Ngcomp,1);
            Ag=zeros(N,Ngcomp);
            bg=obj.p;
            for i=1:Ngcomp
                xx=obj.X(idx==i,:) ;
                ww=obj.p(idx==i);
                [mcp,pcp]=MeanCov(xx,ww/sum(ww));
                Mcomps{i} = mcp(:);
                Pcomps{i} = pcp;
%                 comps{i} = {mcp(:),pcp};
                Ag(:,i)=mvnpdf(obj.X,mcp(:)',pcp);
            end
            
            cvx_begin
                variables wg(Ngcomp) t(N)
                minimize( norm(wg,1)+ErrWt*norm(t,2) )
                subject to
                Ag*wg==bg+t
                wg>0
                %     sum(wg)==1
            cvx_end
            
            GMM.w = wg;
            GMM.mx = Mcomps;
            GMM.Px = Pcomps;
        end
        
    end
end



