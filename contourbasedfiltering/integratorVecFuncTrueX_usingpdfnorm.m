function I = integratorVecFuncTrueX_usingpdfnorm(pdfnorm,VecFunctrue,method)
% integrate Functrue(x) wrt to pdfnorm. Here x is the true space, but given is
% \int Functrue(xtrue)p(xtrue)dxtrue
% ------------------------------------------------------------
% the pdfnorm in norm space
% X are the points that are tracked
% ,LB,UB, GMMhull are all in pdfnorm
dim=pdfnorm.dim;

%estimate normalizing constant


maxdiag = norm(pdfnorm.UB-pdfnorm.LB);
% keyboard
%% method 2 using mixture of gaussians

if strcmp(method,'RegTreeBoxIntegrator')
    %TODO plot the boxes to see if they cover the regions
%     pdfnorm is expected to have pdfnorm.RigTree; which gives the boxes
    boxes=getTree2Boxes(pdfnorm.RigTree.method_params.tree);
    if dim==2
        Nmc=100;
    end
    if dim == 4
    Nmc = 100;
    end
    if dim==6
        Nmc=100;
    end
%     tic
    XX=zeros(Nmc*size(boxes,1),dim);
    VV=zeros(Nmc*size(boxes,1),1);
    k=1;
    for j=1:1:size(boxes,1)
        lb = boxes{j,1};
       ub = boxes{j,2};
       boxdiag = norm(ub-lb);
       volbox = prod(ub-lb );
       
       if boxdiag/maxdiag <= 0.25
           NN=10;
       elseif boxdiag/maxdiag > 0.25 && boxdiag/maxdiag <= 0.5
           NN= 50;
       else
           NN=Nmc;
       end
       xnorm = mvurnd(lb,ub,NN);
       VV(k:k+NN-1)=volbox/NN;
       XX(k:k+NN-1,:)=xnorm;
       k=k+NN;
    end
    XX=XX(1:k-1,:);
    VV=VV(1:k-1);
    I=zeros(length(VecFunctrue),1);
    GG = pdfnorm.func(XX);
    for fn=1:length(VecFunctrue)
        I(fn) = sum(VV.*( VecFunctrue{fn}(pdfnorm.transForms.normX2trueX(XX) ).*GG ));
    end
%     I = sum(VV.*( Functrue(pdfnorm.transForms.normX2trueX(XX) ).*pdfnorm.func(XX) ));
%     toc
    
%     keyboard
    
%     tic
%     I=zeros(size(boxes,1),1);
%     for j=1:1:size(boxes,1)
% %         [j,size(boxes,1)]
%        lb = boxes{j,1};
%        ub = boxes{j,2};
%        boxdiag = norm(ub-lb);
%        volbox = prod(ub-lb );
%        
%        xnorm = mvurnd(lb,ub,Nmc);
%        w=ones(Nmc,1)/Nmc;
%         
%         pdfnorm.func(xnorm);
%         
%        I(j)=sum(w.*( Functrue(pdfnorm.transForms.normX2trueX(xnorm) ).*pdfnorm.func(xnorm) ))*volbox;
%         
%        
%     end
%     I = sum(I);
%     toc
    disp(['Integration constant is :'])
    I
    return 
end




error('Method is not Known/Implemented')






end