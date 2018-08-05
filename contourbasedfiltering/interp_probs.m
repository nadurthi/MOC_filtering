%% max entropic fitting
X0=Xmc;
N=size(Xmc,1);

p0=mvnpdf(X0,[5,5],0.2*eye(2));

X=[Xt(:,61,1),Xt(:,61,2)];
pt=p0;


[m,P]=MeanCov(X,ones(N,1)/N);

Xn=zeros(size(X));
pn=zeros(size(pt));

A=inv(sqrtm(P));
for i=1:N
    Xn(i,:)=A*(X(i,:)-m')';
    
end
pn=pt/det(A)



figure
plot(Xn(:,1),Xn(:,2),'ro')

figure
plot3(Xn(:,1),Xn(:,2),log(pn),'b+')


%% fitting poly to log of probas
dim =2;
Pf=Basis_polyND(dim,4);
% now interpolate the polynomials to Xn and logpn
A=zeros(N,length(Pf));
for r=1:1:N
   A(r,:) = evaluate_MatrixOfPolys(Pf,Xn(r,:));
end
lam = A\log(pn);
mxentpoly=[0,0,0];
for i=1:length(lam)
    mxentpoly=add_sub_polyND(mxentpoly, scalar_multiply_polyND(lam(i),Pf{i}),'add');
end

[x1,x2]=meshgrid(linspace(-2,2,100),linspace(-1,2,100));
pgrid=zeros(size(x1));
for i=1:size(x1,1)
    for j=1:size(x1,2)
        pgrid(i,j) = exp(evaluate_polyND(mxentpoly,[x1(i,j),x2(i,j)]));
    end
end

figure
mesh(x1,x2,pgrid)
hold on
plot3(Xn(:,1),Xn(:,2),pn,'b+','MarkerSize',6)

figure
contour(x1,x2,pgrid)
