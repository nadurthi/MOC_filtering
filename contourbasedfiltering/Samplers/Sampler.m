
classdef Sampler < handle
    % All possible sampling/re-sampleing strategies
    properties
        active_method=0;
        method_params;
        
        possible_methods={'GH',... % do a global polynomial fit tna then take exp
            'GL',... % do decision tree/CART regression
            'GMM-GH',... % do decision tree/CART regression, use GMMHull to add points and adapt the trees
            'GMM-MC',... % just take mean of k-nearest neigbors
            'GMM-GL',... % get nearest neigbors and do smooth poly regression
            'Tree-GL'... % do a cart and then apply gmm to each bin
            };
        
    end
    methods
        % constructor. save a copy constants
        function obj= Sampler(method)
            obj.active_method=method;
        end
        
        function obj= posteriorSamples(obj,Xprior,priorprobs,GMMhull)
            obj.method_params.GMMhull = GMMhull;
            
            if strcmp(obj.active_method, 'GH' )
                obj.method_params.mxentpoly = fitExpPoly_A_Atop_Aineq(X,p,Pf,Xineq,XtestoutsideHull);
            end

            
        end
        
      
        
        
    end
end
