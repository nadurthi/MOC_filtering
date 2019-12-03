function y=evalfuncvec(func,x)
[r,c]=size(x);
if r==1 || c==1
    x=x(:)';
end
y=zeros(size(x,1),1);
for i=1:1:size(x,1)
    y(i) = func(x(i,:));
end