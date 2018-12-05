function plotNDboxes(tree,states)
boxes=getTree2Boxes(tree);
D=[];
for r=1:size(boxes,1)
    mnb = boxes{r,1};
    mxb = boxes{r,2};
    dr=getALLcorners(mnb,mxb);
    D=[D;max(abs(mxb-mnb))];
%         plot3(dsX.X(:,s),dsX.X(:,s+1),pp(:),'bo')
%         hold on
        plot(dr(:,states(1)),dr(:,states(2)),'r','linewidth',2)
%         patch(dr(:,states(1)),dr(:,states(2)),(boxes{r,4}*ones(size(dr,1),1)),'r','linewidth',2)
        axis([-2,2,-2,2])
        title(num2str(r))
%         hold off

    pause(0.2)
    
end