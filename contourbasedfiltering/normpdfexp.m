function y=normpdfexp(x,m,P)
[r,c]=size(x);
if r==1 || c==1
    x=x(:)';
end
invP=inv(P);
y=zeros(size(x,1),1);
c=-log(sqrt(det(2*pi*P)));
for i=1:1:size(x,1)
    a=x(i,:)'-m(:);
    y(i) = c -0.5*a'*invP*a;
end





end