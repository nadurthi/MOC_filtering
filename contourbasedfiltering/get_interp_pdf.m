function [pdf,pdf0Inorm] = get_interp_pdf(X,probs,mquad,Pquad,Nm)
dim =size(X,2);
% 
% [m,P]=MeanCov(X,probs/sum(probs));
c=1;
m=mquad;
P=Pquad/c;


P=diag(diag(P));
% 
N=size(X,1);

Xn=zeros(size(X));

Psqrt=sqrtm(P);
Psqrt_inv=inv(sqrtm(P));
for i=1:N
    Xn(i,:)=Psqrt_inv*(X(i,:)-m')';
end
pn=probs*det(Psqrt);
ind=pn<1e-70;
pn(ind)=1e-70;

logpn = log(pn);
% 
disp('probs nrom')
[det(Psqrt),min(pn),max(pn),max(pn)/min(pn),cond(Pquad)]
ind = sqrt(sum(Xn.^2,2))<6*sqrt(c);
sum(ind)
% Xn=Xn(ind,:);
% pn=pn(ind);
% N=size(Xn,1);

%% -1 to 1 norm

dim =size(X,2);
% 
% [m,P]=MeanCov(X,probs/sum(probs));
mn = min(X,[],1);
Xn=zeros(size(X));
for i=1:N
    Xn(i,:)=X(i,:)-mn;
end
mx = max(Xn,[],1);
A=diag(2./mx);
for i=1:N
    Xn(i,:)=A*Xn(i,:)'-1;
end


pn=probs/det(A);
ind=pn<1e-70;
pn(ind)=1e-70;

logpn = log(pn);
% 
disp('probs nrom')
[det(Psqrt),min(pn),max(pn),max(pn)/min(pn),cond(Pquad)]
ind = sqrt(sum(Xn.^2,2))<6*sqrt(c);
sum(ind)


%% plottinmg
% Xn and pn are the final ones that are fit
figure(33)
plot3(Xn(:,1),Xn(:,2),log(pn),'ro')
% hold on
% plot_nsigellip(mquad(1:2),Pquad(1:2,1:2),2,'g',2)
% hold off
% [IDX,C] = kmeans(Xn, 5);
% hold on
% for i=1:size(C,1)
%     [m,pR]=MeanCov(Xn(IDX==i,:),logpn(IDX==i)/sum(logpn(IDX==i)));
%     plot_nsigellip(m(1:2),1*pR(1:2,1:2),2,'g',2)
% end
% 
% keyboard
%% fitting poly to log of probas
% keyboard
[~,ind]=sort(pn(:),1,'descend');
pn=pn(ind);
Xn=Xn(ind,:);
Pf=Basis_polyND(dim,Nm);
% now interpolate the polynomials to Xn and logpn
% Nleastfit = 3*length(Pf);

A=zeros(N,length(Pf));
tic
for r=1:1:N
   A(r,:) = evaluate_MatrixOfPolys(Pf,Xn(r,:));
end
lam = A\log(pn);
disp('norm fit')
[rank(A),max(abs(A*lam-log(pn))),min(lam),max(lam)]
toc
% evaluate_polyND(Pf{65},Xn(11,:))
% evaluate_polyND_3(Pf{65},Xn(11,:))

Xt=[];
[IDX,C] = kmeans(Xn, 10);
for i=1:size(C,1)
    if length(logpn(IDX==i))>dim*2
        [m,pR]=MeanCov(Xn(IDX==i,:),logpn(IDX==i)/sum(logpn(IDX==i)));
        if all(eig(pR)>0)
            Xt=[Xt;mvnrnd(m,1^2*pR,50)];
        end
    end
end
NMC=size(Xt,1);
% NMC = 3000;
% Xt = [mvnrnd(zeros(dim,1),0.001^2*eye(dim),NMC/3);mvnrnd(zeros(dim,1),0.1^2*eye(dim),NMC/3);mvnrnd(zeros(dim,1),1.5^2*eye(dim),NMC/3)];
At=zeros(NMC,length(Pf));
tic
for r=1:1:NMC
   At(r,:) = evaluate_MatrixOfPolys(Pf,Xt(r,:));
end

%% reguralization points
rad=max(sqrt(sum(Xn.^2,2)));
Xbnd=3*gen_uniform_grid(5,dim);
% Xbnd=1*Xbnd(sqrt(sum(Xbnd.^2,2))>1.5*rad,:);
[size(Xbnd),length(lam)]
removeind=[];
for i=1:size(Xbnd,1)
   Idx = knnsearch(Xn,Xbnd(i,:),'K',10);
   x = Xn(Idx,:);
   mr = mean(x,1);
   r  = max(sqrt(sum((x - repmat(mr,size(x,1),1)).^2)));
   if norm(Xbnd(i,:)-mr)<5*r
       removeind=[removeind,i];
   end
end
size(Xbnd,1)
Xbnd(removeind,:)=[];
figure(33)
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

%%
[~,ind]=sort(pn(:),1,'descend');
ppn=log(pn(ind));
AAn=A(ind,:);
ppn=ppn(:);

% options = optimoptions('lsqlin','Algorithm','interior-point','Display','iter','MaxIterations',200);
% options = optimoptions('fmincon','Algorithm','interior-point','Display','iter','MaxIterations',200,'MaxFunctionEvaluations',10000,'ConstraintTolerance',1e-6);
% options = optimset('Display','iter','MaxFunctionEvaluations',5000);
% lam2 = fmincon(@(lam)norm(A*lam-log(pn)),lam,[],[],AAn,ppn,[],[],[],options);
% lam = lsqlin(A,log(pn),[],[],AAn(1:floor(N/2),:),ppn(1:floor(N/2)),[],[],lam,options);
% lam = lsqlin(A,log(pn),Dineq,-0*ones(M,1),AAn(1:20,:),ppn(1:20),[],[],lam,options);
% AAn(1:floor(N/2),:)    ppn(1:floor(N/2))
% A*lam2 - log(pn)
% AAn*lam2 - ppn

logpn = log(pn);

lamdim=length(lam);
K = -10*ones(size(Dineq,1),1);
KK=K;
DD=Dineq;
AAAn=AAn(1:15,:);
pppn = ppn(1:15);
lenconstr = length(logpn);

% %working good
cvx_begin
    variables t2(15) t(lenconstr) lam2(lamdim)
    minimize( 0.01*norm(lam2,1)+25*norm(t,2)+100*norm(t2,2))
    subject to
    DD*lam2 <= KK  
    A*lam2==logpn+t
    AAAn*lam2==pppn+t2
cvx_end
  
% w = abs(probs);
% w=w/sum(w);

% cvx_begin
%     variables t(lenconstr) lam2(lamdim)
%     minimize( 0.001*norm(lam2,1)+500*norm(w.*t,2))
%     subject to
%     DD*lam2 <= KK  
%     A*lam2==logpn+t
% %     AAAn*lam2==pppn+t2
% cvx_end

%%
% keyboard
lamsol = lam2;
lgpt=At*lamsol;
% lgpt(lgpt>max(logpn))=-10;
figure(33)
plot3(Xn(:,1),Xn(:,2),log(pn),'ro',Xn(:,1),Xn(:,2),A*lamsol,'b+',Xt(:,1),Xt(:,2),lgpt,'gs')
figure(34)
plot3(Xn(:,1),Xn(:,2),pn,'ro',Xn(:,1),Xn(:,2),exp(A*lamsol),'b+',Xt(:,1),Xt(:,2),exp(lgpt),'gs')
% plot_nsigellip(,1,'r',2);




%%
mxentpoly_norm=zeros(1,dim+1);
for i=1:length(lam)
    mxentpoly_norm=add_sub_polyND(mxentpoly_norm, scalar_multiply_polyND(lam(i),Pf{i}),'add');
end
mxentpoly_norm=simplify_polyND(mxentpoly_norm);

% mxentpoly is the 0-I normalized
% sqrtm(P)*x+mu

mxentpoly=linear_transform_poly(mxentpoly_norm,Psqrt_inv,-Psqrt_inv*m(:));
mxentpoly_norm2=linear_transform_poly(mxentpoly,Psqrt,m(:));


% c=1/det(Psqrt);
cexp = -log(det(Psqrt));
c0 = get_coeff_NDpoly(mxentpoly,zeros(1,dim));
mxentpoly = update_or_insert_coeff_NDpoly(mxentpoly,zeros(1,dim),c0+cexp);

pdf.func=@(x)exp(evaluate_polyND(mxentpoly,x));
pdf.poly=mxentpoly;

pdf0Inorm.func=@(x)exp(evaluate_polyND(mxentpoly_norm,x));
pdf0Inorm.poly=mxentpoly_norm;


% keyboard

