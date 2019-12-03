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
for ii=1:size(XMC,1)
    Xmctest2(ii,:) = XMC(ii,:,2);
    Xmctest3(ii,:) = XMC(ii,:,3);
    Xmctest7(ii,:) = XMC(ii,:,9);
end

OE = rv2elm_a(Xmctest2(100,1:6),1)

[N,dim] =size(Xmctest2);

[mXk2,PXk2]=MeanCov(Xmctest2,1/N*ones(N,1));
[mXk3,PXk3]=MeanCov(Xmctest3,1/N*ones(N,1));
[mXk7,PXk7]=MeanCov(Xmctest7,1/N*ones(N,1));

[a,b]=max(pMC);
mXk2=Xmctest2(b,:)';
mXk3=Xmctest3(b,:)';
mXk7=Xmctest7(b,:)';

PXk2=cov(Xmctest2-repmat(mXk2',N,1));
PXk3=cov(Xmctest3-repmat(mXk3',N,1));
PXk7=cov(Xmctest7-repmat(mXk7',N,1));

close all
s=1;
figure
plot(Xmctest2(:,s+1),Xmctest2(:,s+2),'ro')
figure
plot(Xmctest3(:,s+1),Xmctest3(:,s+2),'ro')
figure
plot(Xmctest7(:,s+1),Xmctest7(:,s+2),'ro')

%%
probs = 1/N*ones(N,1);
dsX2 = DataSet(Xmctest2,probs,'TrueState');
dsX3 = DataSet(Xmctest3,probs,'TrueState');
dsX7 = DataSet(Xmctest7,probs,'TrueState');

A=PXk2;
% A(1:3,4:6)=0;
% A(4:6,1:3)=0;

B=PXk3;
% B(1:3,4:6)=0;
% B(4:6,1:3)=0;

C=PXk7;
% C(1:3,4:6)=0;
% C(4:6,1:3)=0;

dsX2.AddMeanCov_to_OI_Trasform(mXk2,3^2*A);
dsX3.AddMeanCov_to_OI_Trasform(mXk3,3^2*B);
dsX7.AddMeanCov_to_OI_Trasform(mXk7,3^2*C);

transForms_at2 = dsX2.GetTrasnformers();
transForms_at3 = dsX3.GetTrasnformers();
transForms_at7 = dsX7.GetTrasnformers();

close all
s=3;
figure
plot(dsX2.X(:,s+1),dsX2.X(:,s+2),'ro')
figure
plot(dsX3.X(:,s+1),dsX3.X(:,s+2),'ro')
figure
plot(dsX7.X(:,s+1),dsX7.X(:,s+2),'ro')



%% 3 to 2
X3=mvurnd(-1*ones(1,dim),1*ones(1,dim),6000);

X3_2=zeros(size(X3));
X3tr=zeros(size(X3));
X2tr=zeros(size(X3));

for i=1:size(X3,1)
    xnorm = X3(i,:);
    xk3=transForms_at3.normX2trueX(xnorm) ;
    X3tr(i,:)=xk3;
    
   xk2 = model.fback(dtkk1,time.Tvec(3),xk3); 
   X2tr(i,:)=xk2;
   
   X3_2(i,:) = transForms_at2.trueX2normX(xk2) ;
   
end
Xmctest2norm = zeros(size(Xmctest2)); 
for i=1:size(Xmctest2,1)
    Xmctest2norm(i,:) =  transForms_at2.trueX2normX(Xmctest2(i,:)) ;
end

close all
for i=1:6
    for j=i+1:6
        figure
        subplot(1,2,1)
        plot(dsX3.X(:,i),dsX3.X(:,j),'ro',X3(:,i),X3(:,j),'b.')
        xlabel(num2str(i));
        ylabel(num2str(j));
        title(['dsX3.X and X3norm'])
        
        subplot(1,2,2)
        plot(dsX2.X(:,i),dsX2.X(:,j),'ro',X3_2(:,i),X3_2(:,j),'b.') %,Xmctest2norm(:,i),Xmctest2norm(:,j),'gs'
        xlabel(num2str(i));
        ylabel(num2str(j));
        title(['dsX2.X and X3_2nomr'])
        
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
%% 3 to 2 from true
% X3=mvurnd(-0.5*ones(1,dim),0.5*ones(1,dim),6000);
X3tr=mvnrnd(mXk3,1^2*PXk3,15000);
X3_2norm=zeros(size(X3tr));
X3norm=zeros(size(X3tr));
X2tr=zeros(size(X3tr));

for i=1:size(X3tr,1)

   xk2tr = model.fback(dtkk1,time.Tvec(3),X3tr(i,:)); 
   X2tr(i,:)=xk2tr;
   
   X3_2norm(i,:) = transForms_at2.trueX2normX(xk2tr) ;
   
   X3norm(i,:) = transForms_at3.trueX2normX(X3tr(i,:)) ;
end


close all
for i=1:6
    for j=i+1:6
        figure
        subplot(1,2,1)
        plot(dsX2.X(:,i),dsX2.X(:,j),'ro',X3_2norm(:,i),X3_2norm(:,j),'b.') %,Xmctest2norm(:,i),Xmctest2norm(:,j),'gs'
        xlabel(num2str(i));
        ylabel(num2str(j));
        title(['dsX2.X and X3_2nomr'])
        subplot(1,2,2)
        plot(dsX3.X(:,i),dsX3.X(:,j),'ro',X3norm(:,i),X3norm(:,j),'b.') %,Xmctest2norm(:,i),Xmctest2norm(:,j),'gs'
        xlabel(num2str(i));
        ylabel(num2str(j));
        title(['dsX3.X and X3norm'])
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
