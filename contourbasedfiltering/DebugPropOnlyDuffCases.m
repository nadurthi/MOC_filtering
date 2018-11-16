clear
close all
clc

load('DuffPropOnlyCasesTest')
savePriorProps.saveit=0;
%%
% test cases 13,16,20,26,27,29,32,34,39
cases = [13,16,20,26,27,29,32,34,39];

k=39;

X=histXprior{k,1};
probs=histXprior{k,2};
    
[mX,PX]=MeanCov(X,probs/sum(probs));

fullnormpdf=get_interp_pdf_0I_duff(X,probs,mX,PX,4,3,k,Xmctest,Xtruth(k,:)); %Xtruth(k,:)
plotpdfs_prior_2D(k,fullnormpdf,X,probs,xfquad,Pfquad,GMM,Xmctest,Xtruth(k,:),savePriorProps)


%%
figure
plot3(X(:,1),X(:,2),probs,'ro')
close all

fullnormpdf=get_interp_pdf_0I_duff(X,probs,mX,PX,4,3,k,Xmctest,Xtruth(k,:)); %Xtruth(k,:)
%     fullnormpdf=get_interp_pdf_0I_2D(X,probs,mX,PX,4,k,[],Xtruth(k,:),plotsconf); %Xtruth(k,:)
priorpdfnorm=fullnormpdf;

plotpdfs_prior_2D(k,fullnormpdf,X,probs,xfquad,Pfquad,GMM,Xmctest,Xtruth(k,:),savePriorProps)

    


