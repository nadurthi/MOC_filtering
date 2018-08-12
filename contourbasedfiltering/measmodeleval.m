function f=measmodeleval(z,x,model)
if any(size(x)==1)
   x=x(:)'; 
end
z=z(:)';
[n,dim] = size(x);
f=zeros(n,1);
% keyboard
for i=1:n
    mm=model.h(x(i,:)');
    f(i) = mvnpdf(z,mm(:)',model.R);
end