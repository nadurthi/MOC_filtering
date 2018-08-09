function [pdfnorm,pdftransF] = get_interp_pdf_hypercube11(X,probs,mquad,Pquad,Nm,Tk,Xmctest)
%% -1 to 1 norm

dim =size(X,2);

N=size(X,1);
% 
% [m,P]=MeanCov(X,probs/sum(probs));
mn = min(X,[],1);
Xn=zeros(size(X));
for i=1:N
    Xn(i,:)=X(i,:)-mn;
end
mx = max(Xn,[],1);
Atransf=diag(2./mx);
mulin = -2*mn(:)./mx(:)-1;
for i=1:N
    Xn(i,:)=Atransf*Xn(i,:)'-1;
end

detAtransf = det(Atransf);

pn=probs/detAtransf;
ind=pn<1e-70;
pn(ind)=1e-70;

logpn = log(pn);



%% plottinmg

figure(33)
plot3(Xn(:,1),Xn(:,2),log(pn),'ro')

%% fitting poly to log of probas
% keyboard
[~,ind]=sort(pn(:),1,'descend');
pn=pn(ind);
Xn=Xn(ind,:);
Pf=Basis_polyND(dim,Nm);
% now interpolate the polynomials to Xn and logpn
% Nleastfit = 3*length(Pf);

A=zeros(N,length(Pf));

for r=1:1:N
   A(r,:) = evaluate_MatrixOfPolys(Pf,Xn(r,:));
end
lam = A\log(pn);
disp('norm fit')
[rank(A),max(abs(A*lam-log(pn))),min(lam),max(lam)]




%% reguralization points
rad=max(sqrt(sum(Xn.^2,2)));
% Xbnd=3*gen_uniform_grid(6,dim);
% Xbnd = 2*(rand(600,dim)*2-1);
% Xbnd=1*Xbnd(sqrt(sum(Xbnd.^2,2))>1.5*rad,:);
[Xbnd,~] = GLgn_pts(-3*ones(1,dim),3*ones(1,dim),7);

[size(Xbnd),length(lam)]
removeind=[];
for i=1:size(Xbnd,1)
   Idx = knnsearch(Xn,Xbnd(i,:),'K',10);
   x = Xn(Idx,:);
   mr = mean(x,1);
   r  = max(sqrt(sum((x - repmat(mr,size(x,1),1)).^2)));
   if norm(Xbnd(i,:)-mr)<2*r
       removeind=[removeind,i];
   end
end
size(Xbnd,1)
Xbnd(removeind,:)=[];
figure(35)
plot3(Xn(:,1),Xn(:,2),log(pn),'ro')
hold on
plot3(Xbnd(:,1),Xbnd(:,2),-ones(size(Xbnd,1),1),'b+')
hold off

% 
M=size(Xbnd,1);
Dineq = zeros(M,length(Pf));
for r=1:1:M
   Dineq(r,:) = evaluate_MatrixOfPolys(Pf,Xbnd(r,:));
end

%% optimization : regularization

[~,ind]=sort(pn(:),1,'descend');
logppn=log(pn(ind));
AAn=A(ind,:);
logppn=logppn(:);


% options = optimoptions('fmincon','Algorithm','interior-point','Display','iter','MaxIterations',200,'MaxFunctionEvaluations',10000,'ConstraintTolerance',1e-6);
% lam2 = fmincon(@(lam)norm(A*lam-log(pn)),lam,[],[],AAn,ppn,[],[],[],options);

logpn = log(pn);

lamdim=length(lam);
K = -5*ones(size(Dineq,1),1);
KK=K;
DD=Dineq;
AAAn=AAn(1:15,:);
logpppn = logppn(1:15);
lenconstr = length(logpn);


% %working good
cvx_begin
    variables t2(15) t(lenconstr) lam2(lamdim)
    minimize( 5*norm(lam2,1)+25*norm(t,2)+100*norm(t2,2))
    subject to
    DD*lam2 <= KK  
    A*lam2==logpn+t
    AAAn*lam2==logpppn+t2
cvx_end
         

lamsol = lam2;

%% normalizing and constructing normalized pdf


mxentpoly_norm=zeros(1,dim+1);
for i=1:length(lamsol)
    mxentpoly_norm=add_sub_polyND(mxentpoly_norm, scalar_multiply_polyND(lamsol(i),Pf{i}),'add');
end
mxentpoly_norm=simplify_polyND(mxentpoly_norm);

pdfnorm.func=@(x)exp(evaluate_polyND(mxentpoly_norm,x));
pdfnorm.polyeval=@(x)evaluate_polyND(mxentpoly_norm,x);
pdfnorm.poly=mxentpoly_norm;
pdfnorm.type = 'true-hypercube-11';


pdfnorm = normalize_exp_pdf(pdfnorm,Xn,mquad,Pquad,'GMM_MC');


pdftransF.trueX2normY = @(x)affineTransform(x,Atransf,mulin(:));
pdftransF.normX2trueX = @(xn)affineTransform(xn,inv(Atransf),inv(Atransf)*mulin(:));
pdftransF.normprob2trueprob = @(p)p/det(inv(Atransf));
pdftransF.trueprob2normprob = @(p)p/detAtransf;
pdftransF.minx = mn;
pdftransF.maxx = mx;
pdftransF.Pquad = Pquad;
pdftransF.mquad = mquad;
pdftransF.Atransf=Atransf;
pdftransF.mulin=mulin;
% pdftransF.Xnorm0I2Xnorm11 = @(xn)Xnorm0I2Xnorm11(xn,minx,maxx,mquad,Pquad);
% pdftransF.Xnorm112Xnorm0I = @(xn)Xnorm112Xnorm0I(xn,minx,maxx,mquad,Pquad);


%% generating test points
Xt=[];
[IDX,C] = kmeans(Xn, 10);
for i=1:size(C,1)
    if length(logpn(IDX==i))>dim*2
        [m,pR]=MeanCov(Xn(IDX==i,:),logpn(IDX==i)/sum(logpn(IDX==i)));
        if all(eig(pR)>0)
            Xt=[Xt;mvnrnd(m,2^2*pR,100)];
        end
    end
end
NMC=size(Xt,1);

At=zeros(NMC,length(Pf));
tic
for r=1:1:NMC
   At(r,:) = evaluate_MatrixOfPolys(Pf,Xt(r,:));
end

lgpt=pdfnorm.polyeval(Xt);
lgpnest=pdfnorm.polyeval(Xn);

figure(33)
plot3(Xn(:,1),Xn(:,2),log(pn),'ro',Xn(:,1),Xn(:,2),lgpnest,'b+',Xt(:,1),Xt(:,2),lgpt,'gs')
title(['time step = ',num2str(Tk)])
figure(34)
plot3(Xn(:,1),Xn(:,2),pn,'ro',Xn(:,1),Xn(:,2),exp(lgpnest),'b+',Xt(:,1),Xt(:,2),exp(lgpt),'gs')
title(['time step = ',num2str(Tk)])
% plot_nsigellip(,1,'r',2);



%% marginal 0I 2D plots
% keyboard

plotmargs=1;

if plotmargs == 1
    [Xx,Xy]=meshgrid(linspace(-1.5,1.5,25),linspace(-1.5,1.5,25) );
    % Xp=[reshape(Xx,625,1),reshape(Xy,625,1)];
    margprobs = zeros(size(Xx));
    for i=1:size(Xx,1)
        for j=1:size(Xx,2)
            margprobs(i,j) = get_2Dmarginalized_probs([Xx(i,j),Xy(i,j)],1,2,Xn,pn,NaN,NaN,pdfnorm,'ClusterMC');
        end
    end


    figure(1)
    contour(Xx,Xy,margprobs,15)
    hold on
    plot(Xn(:,1),Xn(:,2),'ro')
    plot(Xt(:,1),Xt(:,2),'g*')
    title(['time step = ',num2str(Tk)])
    axis equal
    axis square
    hold off
    
    figure(2)
    surf(Xx,Xy,margprobs)
    alpha 0.4
    hold on
    plot(Xn(:,1),Xn(:,2),'ro')
    plot(Xt(:,1),Xt(:,2),'g*')
    title(['time step = ',num2str(Tk)])
    axis equal
    axis square
    hold off
    

end



%% true pdf
% mxentpoly is the 0-I normalized

% 
% mxentpoly=linear_transform_poly(mxentpoly_norm,Atransf,mulin);
% % mxentpoly_norm2=linear_transform_poly(mxentpoly,Psqrt,m(:));
% 
% 
% % c=1/det(Psqrt);
% cexp = log(det(Atransf));
% c0 = get_coeff_NDpoly(mxentpoly,zeros(1,dim));
% mxentpoly = update_or_insert_coeff_NDpoly(mxentpoly,zeros(1,dim),c0+cexp);
% 
% pdf.func=@(x)exp(evaluate_polyND(mxentpoly,x));
% pdf.poly=mxentpoly;
% pdf.type = 'true';

%% final transform functions

% pdftransF.trueX2normX = @(x)Xtrue2Xnorm_OI(x,Pquad,mquad);
% pdftransF.normX2trueX = @(xn)Xnorm2Xtrue_OI(xn,Pquad,mquad);
% pdftransF.normprob2trueprob = @(p)p*detAtransf;
% pdftransF.trueprob2normprob = @(p)p/detAtransf;
% pdftransF.minx = mn;
% pdftransF.maxx = mx;
% pdftransF.Pquad = Pquad;
% pdftransF.mquad = mquad;
% pdftransF.Xnorm0I2Xnorm11 = @(xn)Xnorm0I2Xnorm11(xn,minx,maxx,mquad,Pquad);
% pdftransF.Xnorm112Xnorm0I = @(xn)Xnorm112Xnorm0I(xn,minx,maxx,mquad,Pquad);


