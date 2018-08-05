function [pdf,pdf0Inorm] = get_interp_pdf_nornorm(X,probs,mquad,Pquad,Nm)
dim =size(X,2);
% 
% [m,P]=MeanCov(X,probs/sum(probs));
m=mquad;
P=Pquad;

N=size(X,1);


% Xn and pn are the final ones that are fit

% keyboard
%% fitting poly to log of probas

Pf=Basis_polyND(dim,Nm);
% now interpolate the polynomials to Xn and logpn
A=zeros(N,length(Pf));
for r=1:1:N
   A(r,:) = evaluate_MatrixOfPolys(Pf,X(r,:));
end
lam = A\log(probs);

keyboard

mxentpoly=zeros(1,dim+1);
for i=1:length(lam)
    mxentpoly=add_sub_polyND(mxentpoly, scalar_multiply_polyND(lam(i),Pf{i}),'add');
end
mxentpoly=simplify_polyND(mxentpoly);

pdf.func=@(x)exp(evaluate_polyND(mxentpoly,x));
pdf.poly=mxentpoly;

pdf0Inorm.func=@(x)exp(evaluate_polyND(mxentpoly_norm,x));
pdf0Inorm.poly=mxentpoly_norm;


% keyboard

