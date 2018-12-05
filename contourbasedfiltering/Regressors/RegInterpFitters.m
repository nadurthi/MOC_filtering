
classdef RegInterpFitters < handle
    % All possible fitting/regressor/interpolations
    % Local and global
    properties
        X;
        p;
        poly;
        PfBasis;
        mxentpoly;
        active_method;
        possible_methods={'PolyFit',... % do a global polynomial fit tna then take exp
            'PolyFitLeastSquares', ... % do a least squares fit with l2-norm regularization (analytical)
            'DecisionTree',... % do decision tree/CART regression
            'DecisionTreeAdaptiveOutRegion',... % do decision tree/CART regression, use GMMHull to add points and adapt the trees
            'KnnMean',... % just take mean of k-nearest neigbors
            'KnnMeanKdTree',... % build the kdtree and save it
            'KnnPolyFit',... % get nearest neigbors and do smooth poly regression
            'KnnLeastDeg',... % do k-nn and then least degree interpolation
            'CARTgmm'... % do a cart and then apply gmm to each bin
            };
        method_params;
    end
    methods
        % constructor. save a copy constants
        function obj= RegInterpFitters(method)
            obj.active_method=method;
        end
        
        function obj= fit(obj,X,p,Xineq,XtestoutsideHull,Pf,GMMhull)
            obj.method_params.GMMhull = GMMhull;
            
            if strcmp(obj.active_method, 'PolyFit' )
                obj.method_params.mxentpoly = fitExpPoly_A_Atop_Aineq(X,p,Pf,Xineq,XtestoutsideHull);
            end
            if strcmp(obj.active_method, 'PolyFitLeastSquares' )
                obj.method_params.mxentpoly = fitExpPoly_A_leastsquares(X,p,Pf,Xineq,XtestoutsideHull);
            end

            if strcmp(obj.active_method, 'DecisionTree' )
                obj.method_params.tree = fitTree(X,p,Pf,Xineq,XtestoutsideHull);
            end
            if strcmp(obj.active_method, 'DecisionTreeAdaptiveOutRegion' )
                obj.method_params.tree = fitTree_adaptive2Hull(X,p,Pf,Xineq,XtestoutsideHull,GMMhull);
            end
            if strcmp(obj.active_method, 'CARTgmm' )
                obj.method_params.tree = fitTree(X,p,Pf,Xineq,XtestoutsideHull);
                [GMM,boxes] = fitGMM2tree(X,p,obj.method_params.tree);
                obj.method_params.GMM = GMM;
                obj.method_params.boxes = boxes;
            end
            
            if strcmp(obj.active_method, 'KnnMean' )
                obj.method_params.X = X;
                obj.method_params.Xineq = Xineq;
                obj.method_params.p = p;
            end
             if strcmp(obj.active_method, 'KnnMeanKdTree' )
                [Mdl,Ytrain] = fitkdtree(X,p,Xineq);
                obj.method_params.Mdl = Mdl;
                obj.method_params.Ytrain = Ytrain;
             end
            
            if strcmp(obj.active_method, 'KnnLeastDeg' )
                obj.method_params.X = X;
                obj.method_params.Xineq = Xineq;
                obj.method_params.p = p;
            end
            if strcmp(obj.active_method, 'KnnPolyFit' )
                obj.method_params.X = X;
                obj.method_params.Xineq = Xineq;
                obj.method_params.p = p;
                obj.method_params.Pf = Pf;
            end
            
        end
        %-----------------------------------------------------------------
        function [f,flog]=evalfit(obj,Xtest)
            if strcmp(obj.active_method, 'PolyFit' )
                %                 flog = evaluate_polyND(obj.method_params.mxentpoly,Xtest);
                flog = evaluate_PolyforLargetSetX(obj.method_params.mxentpoly,Xtest);
                f = exp(flog);
            end
             if strcmp(obj.active_method, 'PolyFitLeastSquares' )
                %                 flog = evaluate_polyND(obj.method_params.mxentpoly,Xtest);
                flog = evaluate_PolyforLargetSetX(obj.method_params.mxentpoly,Xtest);
                f = exp(flog);
            end
            if strcmp(obj.active_method, 'DecisionTree' )
                flog = predict(obj.method_params.tree,Xtest);
                f = exp(flog);
            end
            if strcmp(obj.active_method, 'DecisionTreeAdaptiveOutRegion' )
                flog = predict(obj.method_params.tree,Xtest);
                f = exp(flog);
            end
            if strcmp(obj.active_method, 'KnnMean' )
                flog = evalKnnMean(obj.method_params.X,obj.method_params.p,obj.method_params.Xineq,Xtest);
                f = exp(flog);
            end
            if strcmp(obj.active_method, 'KnnMeanKdTree' )
                flog = evalKnnMean_bykdtree(obj.method_params.Mdl,obj.method_params.Ytrain,Xtest);
                f = exp(flog);
            end    
            if strcmp(obj.active_method, 'KnnLeastDeg' )
                flog = evalKnnLeatPoly(obj.method_params.X,obj.method_params.p,obj.method_params.Xineq,Xtest);
                f = exp(flog);
            end
            if strcmp(obj.active_method, 'KnnPolyFit' )
                flog = evalKnnPolyFit(obj.method_params.X,obj.method_params.p,obj.method_params.Xineq,Xtest,obj.method_params.Pf);
                f = exp(flog);
            end
            if strcmp(obj.active_method, 'CARTgmm' )
                f = evalGMM(obj.method_params.GMM,Xtest);
                flog = log(f);
            end
        end
        %-----------------------------------------------------------------
        function f=evalmetric(obj,Xtest,ptest)
            plogtest = log(ptest);
            [~,plogfit]=obj.evalfit(Xtest);
            a=abs(plogtest-plogfit)./plogtest;
            f(1)=max(a);
            f(2)=mean(a);
        end
        %-----------------------------------------------------------------
        function plot(obj,X,p,Xtest,states,LB,UB)
            [f,flog]=obj.evalfit(X);
            ftrue = p;
            flogtrue = log(p);
            
            
            plot3(X(:,states(1)),X(:,states(2)),ftrue,'ro' )
            hold on
            
            if isempty(Xtest)==0
                [ftest,ftestlog]=obj.evalfit(Xtest);
                %                keyboard
                plot3(Xtest(:,states(1)),Xtest(:,states(2)),ftest,'k*' )
            end
            
            
            plot3(X(:,states(1)),X(:,states(2)),f,'bs' )
            for i =1:3
                Xmc = mvurnd(LB,UB,5000);
                [fmc,flogmc]=obj.evalfit(Xmc);
                indmcbad = fmc>max(ftrue);
                fmc(indmcbad);
                plot3(Xmc(indmcbad,states(1)),Xmc(indmcbad,states(2)),fmc(indmcbad),'g+' )
            end
        end
        
        
    end
end

%% Knn least poly

%% Knn poly fit

%% Knn mean

%% RegTree


%% POLYFIT


%-------------------------------------------------------------------


%%
