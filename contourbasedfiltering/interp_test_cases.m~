clc
clear close all

for iii=11:35
    
    load(['case',num2str(iii),'.mat'])
    [iii,k]
    [mX,PX]=MeanCov(Xquad,wquad);
    disp(['cond = ',num2str(cond(PX))])
    fullpdf=get_interpolated_pdf(X,probs,mX,PX,4);
    
    Xtestmc = mvnrnd(mX(:)',1^2*PX,1000);
    pXtest = fullpdf.func(Xtestmc);
    pX = fullpdf.func(X);
    logprobtestmc =log(pXtest);
    logprobX =log(pX);
    figure(10)
    plot3(X(:,1),X(:,2),log(probs),'ro',X(:,1),X(:,2),logprobX,'b+',Xtestmc(:,1),Xtestmc(:,2),logprobtestmc,'gs')
    hold on
    plot_nsigellip(mX(1:2),PX(1:2,1:2),1,'r',2)
    title(['k = ',num2str(k)])
    hold off
    %     saveas(gcf,['nonorm_by4_gh5_k_',num2str(k)],'png')
    figure(11)
    plot3(X(:,1),X(:,2),probs,'ro',X(:,1),X(:,2),pX,'b+',Xtestmc(:,1),Xtestmc(:,2),pXtest,'gs')
    title(['k = ',num2str(k)])
    
    keyboard
    
end