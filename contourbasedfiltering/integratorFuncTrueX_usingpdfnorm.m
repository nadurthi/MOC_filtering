function I = integratorFuncTrueX_usingpdfnorm(pdfnorm,Functrue,method)
% integrate Functrue(x) wrt to pdfnorm. Here x is the true space, but given is
% \int Functrue(xtrue)p(xtrue)dxtrue
% ------------------------------------------------------------
% the pdfnorm in norm space
% X are the points that are tracked
% ,LB,UB, GMMhull are all in pdfnorm
dim=pdfnorm.dim;

%estimate normalizing constant



% keyboard
%% method 2 using mixture of gaussians

if strcmp(method,'RegTreeBoxIntegrator')
    %TODO plot the boxes to see if they cover the regions
%     pdfnorm is expected to have pdfnorm.RigTree; which gives the boxes
    boxes=getTree2Boxes(pdfnorm.RigTree.method_params.tree);
    if dim==2
        Nmc=100;
    end
    if dim==6
        Nmc=100;
    end
    I=0;
    for j=1:1:size(boxes,1)
       lb = boxes{j,1};
       ub = boxes{j,2};
       volbox = prod(ub-lb );
       
       xnorm = mvurnd(lb,ub,Nmc);
       w=ones(Nmc,1)/Nmc;
       I=I+sum(w.*( Functrue(pdfnorm.transForms.normX2trueX(xnorm) ).*pdfnorm.func(xnorm) ))*volbox;
    end
    
    disp(['Integration constant is :',num2str(I)])
    return 
end

if strcmp(method,'GMMhull')

%     pdfnorm is expected to have pdfnorm.GMMhull; which gives the GMM
%     boxes=getTree2Boxes(pdfnorm.RigTree.method_params.tree);
    if dim==2
        Nmc=10000;
    end
    if dim==6
        Nmc=100000;
    end

    
%     C=[];
%     prvstd=100;
%     for i=1:20
        Xmc1 = random(MyGmm2MatlabGMM(pdfnorm.GMMhull),Nmc);
        I = mean( Functrue(pdfnorm.transForms.normX2trueX(Xmc1) ).*( (pdfnorm.func(Xmc1))./evalGMM(GMM,Xmc1) ) ) ;
%         C = [C,c];
%         sdC = std(C);
%         if abs(sdC - prvstd)/prvstd < 0.2
%             disp(['normalization integral converged with #samples = ',num2str(Nmc)])
%             break
%         end
%         prvstd = sdC;
%     end

    
    disp(['Integration constant is :',num2str(I)])
    return 
end



error('Method is not Known/Implemented')






end