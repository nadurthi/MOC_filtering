clc
close all
clear

% gaussian fit 2D
% test UT, GH2, 

% [x,w]=UT_sigmapoints(mu,P,M);
% [xint,wint] = GH_points(mcent,Pcov,Np);

A=[0,1;-5,0 ];
eig(A)
ST=@(t)expm(A*t);
X0=mvnrnd([3,3],[1,0.5;0.5,1],100);
for T=1:0.5:10

for i=1:1:100
Xt(i,:)=ST(T)*X0(i,:)';
end
figure(1)
plot(X0(:,1),X0(:,2),'b+',Xt(:,1),Xt(:,2),'ro')
axis([-10,10,-10,10])
axis square
axis equal
title(num2str(T))
pause(2)
end
T=6.5;
% xdot=Ax;
% sigprior = [2,-1;-1,4];
H=[1,1];
R=1;

sigk = [1,0.8;0.8,1];
mk = [3,3];

priork1sigk = ST(T)*sigk*ST(T)';
priork1mk = ST(T)*mk';

zk1 = H*priork1mk+sqrtm(R)*randn;
Pzk1 = H*priork1sigk*H'+R;
K=priork1sigk*H'*inv(Pzk1);

postk1mk = priork1mk+K*(zk1-H*priork1mk);
postk1sigk = priork1sigk-priork1sigk*H'*inv(Pzk1)*H*priork1sigk;

figure
hold on
plot_nsigellip(mk,sigk,1,'r',2)
plot_nsigellip(priork1mk,priork1sigk,1,'b',2)
plot_nsigellip(postk1mk,postk1sigk,1,'g',2)

pdfkfunc = @(x)mvnpdf(x,mk,sigk);
priorpdfk1func = @(x)mvnpdf(x,priork1mk,priork1sigk);
postpdfk1func = @(x)mvnpdf(x,postk1mk,postk1sigk);


[xprior,wprior] = GH_points(mprior(:),sigprior,3);

%  [xx,yy]= meshgrid(linspace(3-12,3+12,3));
%  xpriorgrid=[xx(:),yy(:)];

[xpriorgrid,w]=GH_points(mprior(:),1*sigprior,3);
xpriorgrid=xpriorgrid+0;
%  xpriorgrid=xprior;
pprior = priorpdffunc(xpriorgrid);
[ap,bp]=MeanCov(xpriorgrid,pprior/sum(pprior));
[aw,bw]=MeanCov(xprior,wprior);
% figure
% plot(1:length(wprior),pprior/sum(pprior),'bo',1:length(wprior),wprior,'r+')
% legend('probs','ghwts')
[ap,bp]
[aw,bw]



Np=size(xpriorgrid,1);

Pf=Basis_polyND(2,2);
numc = length(Pf);
Aeq=zeros(Np,numc);
for ib=1:length(Pf)
    Aeq(:,ib)=evaluate_PolyforLargetSetX(Pf{ib},xpriorgrid);
end
beq = log(pprior)




cvx_begin
variables teq(Np) c(numc) %teqTop(Ntop)
% variable S(2,2) symmetric
minimize( 0.00001*norm(c,1)+1*norm(teq,2) ) %+500*norm(teqTop,2)
subject to
Aeq*c==beq+teq;
% -[2*c(6),c(4);c(4),2*c(5)]>=0;
% S(1,1)==2*c(6);
% S(1,2)==c(4);
% S(2,2)==2*c(5);
% -S>=0;
cvx_end

estsig = inv([2*c(6),c(4);c(4),2*c(5)])
estmu = -estsig*[c(3);c(2)]

