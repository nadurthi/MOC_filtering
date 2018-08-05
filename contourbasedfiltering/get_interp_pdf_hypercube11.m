function [pdf,pdfnorm11,pdfnorm0I,pdftransF] = get_interp_pdf_hypercube11(X,probs,mquad,Pquad,Nm)
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
% 
disp('probs nrom')
[det(Psqrt),min(pn),max(pn),max(pn)/min(pn),cond(Pquad)]
ind = sqrt(sum(Xn.^2,2))<6*sqrt(c);
sum(ind)

keyboard

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

%% optimization : regularization

[~,ind]=sort(pn(:),1,'descend');
ppn=log(pn(ind));
AAn=A(ind,:);
ppn=ppn(:);


% options = optimoptions('fmincon','Algorithm','interior-point','Display','iter','MaxIterations',200,'MaxFunctionEvaluations',10000,'ConstraintTolerance',1e-6);
% lam2 = fmincon(@(lam)norm(A*lam-log(pn)),lam,[],[],AAn,ppn,[],[],[],options);

logpn = log(pn);

lamdim=length(lam);
K = -10*ones(size(Dineq,1),1);
KK=K;
DD=Dineq;
lenconstr = length(logpn);

cvx_begin
    variables t(lenconstr) lam2(lamdim)
    minimize( 1*norm(lam2,1)+500*norm(t,2))
    subject to
    DD*lam2 <= KK  
    A*lam2==logpn+t
cvx_end
         

%% generating test points
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

At=zeros(NMC,length(Pf));
tic
for r=1:1:NMC
   At(r,:) = evaluate_MatrixOfPolys(Pf,Xt(r,:));
end

lamsol = lam2;
lgpt=At*lamsol;
% lgpt(lgpt>max(logpn))=-10;
figure(33)
plot3(Xn(:,1),Xn(:,2),log(pn),'ro',Xn(:,1),Xn(:,2),A*lamsol,'b+',Xt(:,1),Xt(:,2),lgpt,'gs')
figure(34)
plot3(Xn(:,1),Xn(:,2),pn,'ro',Xn(:,1),Xn(:,2),exp(A*lamsol),'b+',Xt(:,1),Xt(:,2),exp(lgpt),'gs')
% plot_nsigellip(,1,'r',2);




%% norm pdf -- also it should be normalized
mxentpoly_norm=zeros(1,dim+1);
for i=1:length(lam)
    mxentpoly_norm=add_sub_polyND(mxentpoly_norm, scalar_multiply_polyND(lam(i),Pf{i}),'add');
end
mxentpoly_norm=simplify_polyND(mxentpoly_norm);

pdfnorm11.func=@(x)exp(evaluate_polyND(mxentpoly_norm,x));
pdfnorm11.poly=mxentpoly_norm;
pdfnorm11.type = 'hypercube-11';

% normaliziing the norm pdf
pdfnorm11 = normalize_exp_pdf(pdfnorm11,Xn,mquad,Pquad,'GMM_MC');


% getting the 0I norm pdf
P0I=linear_transform_poly(pdfnorm11.poly,Atransf*sqrtm(Pquad),A*mquad(:)+mulin);
pdfnorm0I.func=@(x)exp(evaluate_polyND(P0I,x));
pdfnorm0I.poly=P0I;
pdfnorm0I.type = 'OI';


% marginalize the pdf

%% true pdf
% mxentpoly is the 0-I normalized


mxentpoly=linear_transform_poly(mxentpoly_norm,Atransf,mulin);
% mxentpoly_norm2=linear_transform_poly(mxentpoly,Psqrt,m(:));


% c=1/det(Psqrt);
cexp = log(det(Atransf));
c0 = get_coeff_NDpoly(mxentpoly,zeros(1,dim));
mxentpoly = update_or_insert_coeff_NDpoly(mxentpoly,zeros(1,dim),c0+cexp);

pdf.func=@(x)exp(evaluate_polyND(mxentpoly,x));
pdf.poly=mxentpoly;
pdf.type = 'true';

%% final transform functions

pdftransF.trueX2normX = @(x)Xtrue2Xnorm_OI(x,Pquad,mquad);
pdftransF.normX2trueX = @(xn)Xnorm2Xtrue_OI(xn,Pquad,mquad);
pdftransF.normprob2trueprob = @(p)p*detAtransf;
pdftransF.trueprob2normprob = @(p)p/detAtransf;
pdftransF.minx = mn;
pdftransF.maxx = mx;
pdftransF.Pquad = Pquad;
pdftransF.mquad = mquad;
pdftransF.Xnorm0I2Xnorm11 = @(xn)Xnorm0I2Xnorm11(xn,minx,maxx,mquad,Pquad);
pdftransF.Xnorm112Xnorm0I = @(xn)Xnorm112Xnorm0I(xn,minx,maxx,mquad,Pquad);


