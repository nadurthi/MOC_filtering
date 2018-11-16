function TreeHullGL_Sampler()
tree = fitrtree(XmcHull,Y,'MinParentSize',50,'MaxNumSplits',500,'MinLeafSize',50,...
                    'PredictorNames',{'x','y'},'ResponseName','probs');
boxes=getTree2Boxes(tree);
figure(105)
plot3(XmcHull(:,1),XmcHull(:,2),Y,'ro')
hold on
plot2Dboxes(tree)
hold off
