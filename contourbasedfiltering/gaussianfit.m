

% gaussian fit 2D
% test UT, GH2, 

[x,w]=UT_sigmapoints(mu,P,M);
[xint,wint] = GH_points(mcent,Pcov,Np);

A=[0,1;-5,0 ];
eig(A)
ST=@(t)expm(A*t);
X0=mvnrnd([3,3],0.5*eye(2),100);
for T=1:0.5:10

for i=1:1:100
Xt(i,:)=ST(T)*X0(i,:)';
end
figure(1)
plot(X0(:,1),X0(:,2),'b+',Xt(:,1),Xt(:,2),'ro')

pause(2)
end
% xdot=Ax;
% sigprior = [2,-1;-1,4];
sigprior = [1,0.8;0.8,1];
mprior = [0,0];
priorpdffunc = @(x)mvnpdf(x,mprior,sigprior);

[xprior,wprior] = GH_points(mprior(:),sigprior,3);

%  [xx,yy]= meshgrid(linspace(3-12,3+12,3));
%  xpriorgrid=[xx(:),yy(:)];

% [xpriorgrid,w]=GH_points(mprior(:),1.5*sigprior,3);
xpriorgrid=mvnrnd(mprior(:),1.5*sigprior,15);
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

figure
plot3(xpriorgrid(:,1),xpriorgrid(:,2),pprior,'ro')
hold on
plot_nsigellip(mprior,sigprior,1,'r',2)

Np=size(xpriorgrid,1);

Pf=Basis_polyND(2,2);
numc = length(Pf);
Aeq=zeros(Np,numc);
for ib=1:length(Pf)
    Aeq(:,ib)=evaluate_PolyforLargetSetX(Pf{ib},xpriorgrid);
end
pq=zeros(Np,1);
pq(randi(Np,1))=1;
pq(randi(Np,1))=1;
pq(randi(Np,1))=1;
beq = log(pprior)+1;

err=beq-log(pprior)


cvx_begin
variables teq(Np) c(numc) %teqTop(Ntop)
% variable S(2,2) symmetric
minimize( 0.000001*norm(c,1)+1*norm(teq,2) ) %+500*norm(teqTop,2)
subject to
Aeq*c==beq+teq;
% c(6)<=0;
% c(5)<=0;

% [2*c(6),c(4);c(4),2*c(5)]<=0;
% S(1,1)==2*c(6);
% S(1,2)==c(4);
% S(2,2)==2*c(5);
% -S>=0;
cvx_end

estsig = inv([2*c(6),c(4);c(4),2*c(5)])
estmu = -estsig*[c(3);c(2)]

