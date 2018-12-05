function Xineq = getXineq(RigTree,GMMHull,MaxN,LB,UB)

mainboxdiag = norm(UB-LB);
minlnprob = min(RigTree.method_params.tree.NodeMean);
maxlnprob = max(RigTree.method_params.tree.NodeMean);
lambda = 0.7;
cutofflnprobs = lambda*minlnprob + (1-lambda)*maxlnprob;

boxes=getTree2Boxes(RigTree.method_params.tree,LB,UB);
Xineq=[];
boxfrac_vec=zeros(size(boxes,1),1);

for j=1:size(boxes,1)
    lb = boxes{j,1};
    ub = boxes{j,2};
    
    boxdiag = norm(ub-lb);
    
    
    if boxes{j,4} > cutofflnprobs
        continue
    end
    
    boxfrac = boxdiag/mainboxdiag;
    boxfrac_vec(j) = boxfrac;
    
    if boxfrac<0.4
        s=0.5*(lb+ub);
        Xineq = vertcat(Xineq,s(:)');
        
        continue
    end
    
    if boxfrac>=0.4 && boxfrac<0.6
%         [Xbnd1,~] = GLgn_pts(lb(:)',ub(:)',2);
        Xbnd1=mvurnd(lb,ub,10);
        Xineq = vertcat(Xineq,Xbnd1);
        continue
    end
    
    
    if boxfrac>=0.6 && boxfrac<0.8
%         [Xbnd1,~] = GLgn_pts(lb(:)',ub(:)',3);
        Xbnd1=mvurnd(lb,ub,100);
        Xineq = vertcat(Xineq,Xbnd1);
        continue
    end
    
    if boxfrac>=0.8
%         [Xbnd1,~] = GLgn_pts(lb(:)',ub(:)',4);
        Xbnd1=mvurnd(lb,ub,1000);
        Xineq = vertcat(Xineq,Xbnd1);
        continue
    end
    
end
indbnd = GMMHull.IsInsideHull(Xineq,1.7);

Xineq = Xineq(~indbnd,:);

% keyboard

figure(312)
hist(boxfrac_vec(boxfrac_vec>0),linspace(0,1,100))


