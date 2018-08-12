function f=gausspdfmodeleval(z,x,R)
if any(size(x)==1)
   x=x(:)'; 
end
z=z(:)';
[n,dim] = size(x);
f=zeros(n,1);
for i=1:n
    f(i) = mvnpdf(z,x(i,:),R);
end