function Y=affineTransform(X,A,MU)

if any(size(X)==1)
   X=X(:)'; 
end
dim=size(X,2);
n = size(X,1);
MU=MU(:);

Y=zeros(size(X));
for i=1:n
   Y(i,:) = A*X(i,:)'+MU; 
end
    