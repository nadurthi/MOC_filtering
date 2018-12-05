function boxes=getTree2Boxes(tree,LB,UB)
states=tree.PredictorNames;
dim = length(states);

switch nargin
    case 3
    
    case 2
        error('Not enough arguments')
    case 1
        LB=-1.5*ones(1,dim);
        UB=1.5*ones(1,dim);
    otherwise
        LB=-1.5*ones(1,dim);
        UB=1.5*ones(1,dim);
end
    

boxes={LB,UB,1,0};


for nc = 1:length(tree.CutPredictor)
    st = tree.CutPredictor{nc};
    if strcmp(tree.CutPredictor{nc},'')==false % node has a branch
        cp = tree.CutPoint(nc);
        stind = getindex2cell(states,st);
        for r=1:nc+10
            if boxes{r,3}==nc
                break
            end
        end
        
        boxmin = boxes{r,1};
        boxmax = boxes{r,2};
        
        newminL=boxmin;
        newmaxL = boxmax;
        newmaxL(stind) = cp;
        
        newminR=boxmin;
        newmaxR = boxmax;
        newminR(stind) = cp;
        
        l = size(boxes,1);
        boxes{l+1,1}=newminL;
        boxes{l+1,2}=newmaxL;
        boxes{l+1,3}=tree.Children(nc,1);
        boxes{l+1,4} = tree.NodeMean(boxes{l+1,3});
        
        boxes{l+2,1}=newminR;
        boxes{l+2,2}=newmaxR;
        boxes{l+2,3}=tree.Children(nc,2);
        boxes{l+2,4} = tree.NodeMean(boxes{l+2,3});
        
        boxes(r,:)=[];
        
    else  % node has no brachs or cildern ... hence it is a leaf
        
    end
    
end
