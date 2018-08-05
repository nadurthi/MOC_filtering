function X=Xnorm2Xtrue_hypercube_11(Xn,minx,maxx)
% hypercube is -1 to 1
if any(size(Xn)==1)
   Xn=Xn(:)'; 
end
dim=size(Xn,2);
n = size(Xn,1);

minx = minx(:)';
maxx = maxx(:)';

A=diag(maxx/2);
X = (Xn + 1*ones(n,dim))/2;
X = X.*repmat(maxx,n,1);
X = X + repmat(minx,n,1);







