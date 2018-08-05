function [pdf,pdf0Inorm] = get_interp_pdf_nonlin(X,probs,mquad,Pquad,Nm)
dim =size(X,2);
% 
% [m,P]=MeanCov(X,probs/sum(probs));
m=mquad;
P=Pquad;

N=size(X,1);

Xn=zeros(size(X));

Psqrt=sqrtm(P);
Psqrt_inv=inv(sqrtm(P));
for i=1:N
    Xn(i,:)=Psqrt_inv*(X(i,:)-m')';
end
pn=probs/det(Psqrt_inv);
ind=pn<1e-30;
pn(ind)=1e-30;

figure
plot3(Xn(:,1),Xn(:,2),pn,'ro')
keyboard

% Xn and pn are the final ones that are fit

% keyboard
%% fitting poly to log of probas

Pf=Basis_polyND(dim,Nm);
% now interpolate the polynomials to Xn and logpn
A=zeros(N,length(Pf));
for r=1:1:N
   A(r,:) = evaluate_MatrixOfPolys(Pf,Xn(r,:));
end
lam = A\log(pn);

% keyboard

r=max(sqrt(sum(Xn.^2,2)));
% while(1)
%     Xbnd=mvnrnd(zeros(1,dim),(3*r)^2*eye(dim),5000);
%     Xbnd=Xbnd(sqrt(sum(Xbnd.^2,2))>4*r,:);
%     if size(Xbnd,1)>100
%         break
%     end
% end
% Xbnd = sphere4Dm(5);
Xbnd=[1.2*r*sphere4Dm(5);5*r*sphere4Dm(6)];
Xbnd=6*gen_uniform_grid(10,4);
Xbnd=Xbnd(sqrt(sum(Xbnd.^2,2))>10,:);
% keyboard

% Xbnd=Xbnd(1:1000,:);
M=size(Xbnd,1);
Dineq = zeros(M,length(Pf));
for r=1:1:M
   Dineq(r,:) = evaluate_MatrixOfPolys(Pf,Xbnd(r,:));
end
% keyboard
[~,ind]=sort(pn(:),1,'descend');
ppn=log(pn(ind));
AAn=A(ind,:);
ppn=ppn(:);
options = optimoptions('lsqlin','Algorithm','interior-point','Display','iter','MaxIterations',200);
lam = lsqlin(A,log(pn),Dineq,-10*ones(M,1),AAn(1:floor(N/100),:),ppn(1:floor(N/100)),[],[],lam,options);

% [~,ind]=sort(pn);
% A=A(ind,:);
% ppn=pn(ind);
% W=diag([ones(1,N-floor(N/2)),5*ones(1,floor(N/2))]);
% lam = (A'*W*A)\(A'*W*log(ppn));
% keyboard

mxentpoly_norm=zeros(1,dim+1);
for i=1:length(lam)
    mxentpoly_norm=add_sub_polyND(mxentpoly_norm, scalar_multiply_polyND(lam(i),Pf{i}),'add');
end
mxentpoly_norm=simplify_polyND(mxentpoly_norm);

% mxentpoly is the 0-I normalized
% sqrtm(P)*x+mu

mxentpoly=linear_transform_poly(mxentpoly_norm,Psqrt_inv,-Psqrt_inv*m(:));

% c=1/det(Psqrt);
cexp = -log(det(Psqrt));
c0 = get_coeff_NDpoly(mxentpoly,zeros(1,dim));
mxentpoly = update_or_insert_coeff_NDpoly(mxentpoly,zeros(1,dim),c0+cexp);

pdf.func=@(x)exp(evaluate_polyND(mxentpoly,x));
pdf.poly=mxentpoly;

pdf0Inorm.func=@(x)exp(evaluate_polyND(mxentpoly_norm,x));
pdf0Inorm.poly=mxentpoly_norm;


% keyboard

