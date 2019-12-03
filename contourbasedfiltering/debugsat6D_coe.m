clc
clear all
close all

load('6Dtestcase_backprop.mat')
dim=6;
dtkk1 = time.dt;


kk=3;
Tstepk1=kk;
Tk1=time.Tvec(kk);

% [X,probs]=propagate_character(X,probs,time.dt,time.Tvec(kk),model);
% [aa,bb]=MeanCov(X,probs/sum(probs));

Xmctest2 = zeros(size(XMC,1),model.fn);
Xmctest3 = zeros(size(XMC,1),model.fn);
Xmctest7 = zeros(size(XMC,1),model.fn);

Xmctest2coe = zeros(size(XMC,1),model.fn);
Xmctest3coe = zeros(size(XMC,1),model.fn);
Xmctest7coe = zeros(size(XMC,1),model.fn);

for ii=1:size(XMC,1)
    Xmctest2(ii,:) = XMC(ii,:,2);
    Xmctest3(ii,:) = XMC(ii,:,3);
    Xmctest7(ii,:) = XMC(ii,:,7);
    
    Xmctest2coe(ii,:) = rv2elm_a(Xmctest2(ii,:),1);
    Xmctest3coe(ii,:) = rv2elm_a(Xmctest3(ii,:),1);
    Xmctest7coe(ii,:) = rv2elm_a(Xmctest7(ii,:),1);
    
end

% test reconvert
X2rv=zeros(size(Xmctest2coe));
for ii=1:size(XMC,1)
    X2rv(ii,:) = elm2rv_a(Xmctest2coe(ii,:),0,1);
end

[N,dim] =size(Xmctest2);
[mXk2,PXk2]=MeanCov(Xmctest2,1/N*ones(N,1));
[mXk3,PXk3]=MeanCov(Xmctest3,1/N*ones(N,1));
[mXk7,PXk7]=MeanCov(Xmctest7,1/N*ones(N,1));


[mXk2coe,PXk2coe]=MeanCov(Xmctest2coe,1/N*ones(N,1));
[mXk3coe,PXk3coe]=MeanCov(Xmctest3coe,1/N*ones(N,1));
[mXk7coe,PXk7coe]=MeanCov(Xmctest7coe,1/N*ones(N,1));



close all
s=4;
figure
plot(Xmctest2coe(:,s+1),Xmctest2coe(:,s+2),'ro')
figure
plot(Xmctest3coe(:,s+1),Xmctest3coe(:,s+2),'ro')
figure
plot(Xmctest7coe(:,s+1),Xmctest7coe(:,s+2),'ro')


for i=1:6
    for j=i+1:6
        figure
        plot(Xmctest2coe(:,i),Xmctest2coe(:,j),'ro')
        xlabel(num2str(i));
        ylabel(num2str(j));
    end
end

%%
probs = 1/N*ones(N,1);
dsX2coe = DataSet(Xmctest2coe,probs,'TrueState');
dsX3coe = DataSet(Xmctest3coe,probs,'TrueState');
dsX7coe = DataSet(Xmctest7coe,probs,'TrueState');


dsX2coe.AddMeanCov_to_OI_Trasform(mXk2coe,4^2*PXk2coe);
dsX3coe.AddMeanCov_to_OI_Trasform(mXk3coe,4^2*PXk3coe);
dsX7coe.AddMeanCov_to_OI_Trasform(mXk7coe,4^2*PXk7coe);

transForms_at2coe = dsX2coe.GetTrasnformers();
transForms_at3coe = dsX3coe.GetTrasnformers();
transForms_at7coe = dsX7coe.GetTrasnformers();

close all
s=3;
figure
plot(dsX2coe.X(:,s+1),dsX2coe.X(:,s+2),'ro')
figure
plot(dsX3coe.X(:,s+1),dsX3coe.X(:,s+2),'ro')
figure
plot(dsX7coe.X(:,s+1),dsX7coe.X(:,s+2),'ro')


for i=1:6
    for j=i+1:6
        figure
        plot(dsX2coe.X(:,i),dsX2coe.X(:,j),'ro')
        xlabel(num2str(i));
        ylabel(num2str(j));
    end
end

%% 3 to 2 from true
X3normcoe=mvurnd(-1*ones(1,dim),1*ones(1,dim),6000);
% X3tr=mvnrnd(mXk3,1^2*PXk3,15000);
X3trcoe=zeros(size(X3normcoe));
X3tr=zeros(size(X3normcoe));
X2tr=zeros(size(X3normcoe));
X2trcoe=zeros(size(X3normcoe));
X2normcoe=zeros(size(X3normcoe));


for i=1:size(X3normcoe,1)
    x3trcoe=transForms_at3coe.normX2trueX(X3normcoe(i,:));
    X3tr(i,:) = elm2rv_a(x3trcoe,0,1);
    
    X2tr(i,:) = model.fback(dtkk1,time.Tvec(3),X3tr(i,:));
    
    X2trcoe(i,:) = rv2elm_a(X2tr(i,:),1);
    
    X2normcoe(i,:) = transForms_at2coe.trueX2normX(X2trcoe(i,:)) ;

end


close all
for i=1:6
    for j=i+1:6
        figure
        subplot(1,2,1)
        plot(dsX2coe.X(:,i),dsX2coe.X(:,j),'ro',X2normcoe(:,i),X2normcoe(:,j),'b.') %,Xmctest2norm(:,i),Xmctest2norm(:,j),'gs'
        xlabel(num2str(i));
        ylabel(num2str(j));
        title(['dsX2coe.X and X2normcoe'])
        subplot(1,2,2)
        plot(dsX3coe.X(:,i),dsX3coe.X(:,j),'ro',X3normcoe(:,i),X3normcoe(:,j),'b.') %,Xmctest2norm(:,i),Xmctest2norm(:,j),'gs'
        xlabel(num2str(i));
        ylabel(num2str(j));
        title(['dsX3coe.X and X3normcoe'])
    end
end


close all
for i=1:6
    for j=i+1:6
        figure
        subplot(1,2,1)
        plot(Xmctest3(:,i),Xmctest3(:,j),'ro',X3tr(:,i),X3tr(:,j),'bs')
        xlabel(num2str(i));
        ylabel(num2str(j));
        title(['Xmctest3 and X3tr'])
        
        subplot(1,2,2)
        plot(Xmctest2(:,i),Xmctest2(:,j),'ro',X2tr(:,i),X2tr(:,j),'bs')
        xlabel(num2str(i));
        ylabel(num2str(j));
        title(['Xmctest2 and X2tr'])
        
    end
end
