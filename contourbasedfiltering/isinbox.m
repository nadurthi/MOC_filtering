function insidept = isinbox(x,LB,UB)
[r,c]=size(x);
if r==1 || c==1
    x=x(:)';
end
insidept=zeros(size(x,1),1);
for i=1:size(x,1)
    if all( (x(i,:)'-LB(:))>=0 ) && all( (UB(:)-x(i,:)')>=0 )
        insidept(i)=1;
    end
end