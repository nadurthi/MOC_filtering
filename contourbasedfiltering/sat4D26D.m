function y=sat4D26D(x)
[r,c]=size(x);
if r==1 || c==1
    x=x(:)';
end
[N,dim]=size(x);

y=zeros(N,6);
y(:,[1,2,4,5])=x;

if r==1 || c==1
    y=y(:);
end