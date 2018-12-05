function plot2Dboxes(tree,LB,UB)
boxes=getTree2Boxes(tree,LB,UB);
D=[];
for r=1:size(boxes,1)
    mnb = boxes{r,1};
    mxb = boxes{r,2};
    dr=getALLcorners(mnb,mxb);
    dr=[dr;dr(1,:)];
    D=[D;max(abs(mxb-mnb))];
%         plot3(dsX.X(:,s),dsX.X(:,s+1),pp(:),'bo')
        hold on
%         patch(dr(:,1),dr(:,2),(boxes{r,4}*ones(size(dr,1),1)),'r','linewidth',2)
%         patch(dr(:,1),dr(:,2),exp(boxes{r,4}*ones(size(dr,1),1)),'r','linewidth',2)
            plot(dr(:,1),dr(:,2),'r','linewidth',2)
        axis([-2,2,-2,2])
        title(num2str(r))
%         hold off

    pause(0.2)
    
end