function X=Xnorm2Xtrue_OI(Xn,P,Mu)
% hypercube is -1 to 1
if any(size(Xn)==1)
   Xn=Xn(:)'; 
end
dim=size(Xn,2);
n = size(Xn,1);
Mu=Mu(:);

A=sqrtm(P);

X=zeros(size(Xn));
for i=1:n
    X(i,:) = A*Xn(i,:)'+Mu;
end
