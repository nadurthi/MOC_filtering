function Xnoi=Xnorm112Xnorm0I(Xn11,minx,maxx,mquad,Pquad)
% hypercube is -1 to 1
if any(size(Xn11)==1)
   Xn11=Xn11(:)'; 
end
dim=size(Xn11,2);
n = size(Xn11,1);
mquad=mquad(:);
minx = minx(:)';
maxx = maxx(:)';

A=diag(2./maxx);
MU=-2*minx(:)./maxx(:) -1;

Amultinv = inv(A*sqrtm(Pquad));
Xnoi=zeros(size(Xn11));
for i=1:n
    Xnoi(i,:) = Amultinv*(Xn11(i,:)'-A*mquad-MU);
end








