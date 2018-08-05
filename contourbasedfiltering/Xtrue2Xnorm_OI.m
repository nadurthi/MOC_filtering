function Xn=Xtrue2Xnorm_OI(X,P,Mu)
% hypercube is -1 to 1
if any(size(X)==1)
   X=X(:)'; 
end
dim=size(X,2);
n = size(X,1);
Mu=Mu(:);

A=sqrtm(inv(P));

Xn=zeros(size(X));
for i=1:n
    Xn(i,:) = A*(X(i,:)'-Mu);
end

