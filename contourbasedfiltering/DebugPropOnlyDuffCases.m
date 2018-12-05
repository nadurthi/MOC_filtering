clear
close all
clc

load('DuffPropOnlyCasesTest')
savePriorProps.saveit=0;
%%
% test cases 13,16,20,26,27,29,32,34,39
cases = [13,16,20,26,27,29,32,34,39];

k=6;

X=histXprior{k,1};
probs=histXprior{k,2};
    
[mX,PX]=MeanCov(X,probs/sum(probs));

fullnormpdf=get_interp_pdf_0I_duff(X,probs,mX,PX,4,3,k,Xmctest,Xtruth(k,:)); %Xtruth(k,:)
plotpdfs_prior_2D(k,fullnormpdf,X,probs,xfquad,Pfquad,GMM,Xmctest,Xtruth(k,:),savePriorProps)





%%
tree=fullnormpdf.RigTree.method_params.tree;
LB=fullnormpdf.LB;
UB=fullnormpdf.UB;
boxes=getTree2Boxes(tree,LB,UB)

figure
plot2Dboxes(tree,LB,UB)
xlabel('x_1')
ylabel('x_2')
axis square
axis equal
axis([-2,2,-2,2])
plot_prop_paper44

figure
plot2Dpatches(tree,LB,UB)
xlabel('x_1')
ylabel('x_2')
axis square
axis equal
axis([-2,2,-2,2])
plot_prop_paper44