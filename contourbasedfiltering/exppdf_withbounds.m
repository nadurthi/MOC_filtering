function pdf = exppdf_withbounds(x,maxentpoly,LB,UB)
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
pdf=zeros(size(x,1),1);
pdf(insidept==1) = exp(evaluate_polyND(maxentpoly,x(insidept==1,:)));

% pdf= exp(evaluate_polyND(maxentpoly,x));