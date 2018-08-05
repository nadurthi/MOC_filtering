function Xn=Xtrue2Xnorm_hypercube_11(X,minx,maxx)
% hypercube is -1 to 1
if any(size(X)==1)
   X=X(:)'; 
end
dim=size(X,2);
n = size(X,1);

minx = minx(:)';
maxx = maxx(:)';

A=diag(2./maxx);

Xn = X - repmat(minx,n,1);
Xn = 2*Xn./repmat(maxx,n,1) -1*ones(n,dim);







