function Xn11=Xnorm0I2Xnorm11(Xnoi,minx,maxx,mquad,Pquad)
% hypercube is -1 to 1
if any(size(Xnoi)==1)
   Xnoi=Xnoi(:)'; 
end
dim=size(Xnoi,2);
n = size(Xnoi,1);
mquad=mquad(:);
minx = minx(:)';
maxx = maxx(:)';

A=diag(2./maxx);
MU=-2*minx(:)./maxx(:) -1;

Amult = A*sqrtm(Pquad);
Xn11=zeros(size(Xnoi));
for i=1:n
    Xn11(i,:) = Amult*Xnoi(i,:)'+A*mquad+MU;
end








