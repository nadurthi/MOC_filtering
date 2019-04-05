for k=2:time.Ntsteps
%     close all
    
    
    X=histXprior{k,1};
    probs=histXprior{k,2};
    priorfullnormpdf=histXprior{k,3};
    
    xfquad=histXprior{k,4};
    Pfquad=histXprior{k,5};
    
    plotpdfs_prior_6D([1,2],k,priorfullnormpdf,X,probs,xfquad,Pfquad,Xmctest,Xtruth(k,:),model,savePriorProps)
    plotpdfs_prior_6D([2,3],k,priorfullnormpdf,X,probs,xfquad,Pfquad,Xmctest,Xtruth(k,:),model,savePriorProps)
    plotpdfs_prior_6D([3,4],k,priorfullnormpdf,X,probs,xfquad,Pfquad,Xmctest,Xtruth(k,:),model,savePriorProps)
    plotpdfs_prior_6D([4,5],k,priorfullnormpdf,X,probs,xfquad,Pfquad,Xmctest,Xtruth(k,:),model,savePriorProps)
    plotpdfs_prior_6D([5,6],k,priorfullnormpdf,X,probs,xfquad,Pfquad,Xmctest,Xtruth(k,:),model,savePriorProps)
    
%     keyboard
    
    
    
%    keyboard
    pause(1)
    
    
    

             X=histXpost{k,1};
            probs=histXpost{k,2};
            fullnormpdf=histXpost{k,3};

            xfquad=histXpost{k,4};
            Pfquad=histXpost{k,5};
            zk=histXpost{k,6};
            plotpdfs_post_6D([1,2],k,priorfullnormpdf,fullnormpdf,model,X,probs,xfquad,Pfquad,Xmctest,Xtruth(k,:),zk,savePostProps)
            plotpdfs_post_6D([2,3],k,priorfullnormpdf,fullnormpdf,model,X,probs,xfquad,Pfquad,Xmctest,Xtruth(k,:),zk,savePostProps)
            plotpdfs_post_6D([3,4],k,priorfullnormpdf,fullnormpdf,model,X,probs,xfquad,Pfquad,Xmctest,Xtruth(k,:),zk,savePostProps)
            plotpdfs_post_6D([4,5],k,priorfullnormpdf,fullnormpdf,model,X,probs,xfquad,Pfquad,Xmctest,Xtruth(k,:),zk,savePostProps)
            plotpdfs_post_6D([5,6],k,priorfullnormpdf,fullnormpdf,model,X,probs,xfquad,Pfquad,Xmctest,Xtruth(k,:),zk,savePostProps)
            


end