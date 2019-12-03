function pdfexp = expofaexppdf_withbounds(x,maxentpoly,LB,UB)


[r,c]=size(x);


if r==1 || c==1
    x=x(:)';
end

[N,dim]=size(x);
LB=-1.0*ones(1,dim);
UB=1.0*ones(1,dim);

insidept=zeros(size(x,1),1);
for i=1:size(x,1)
    if all( (x(i,:)'-LB(:))>=0 ) && all( (UB(:)-x(i,:)')>=0 )
        insidept(i)=1;
    end
end
% m = 0.5*(UB+LB);
% P = diag((0.01*(UB-LB)).^2);
M=-50;
pdfexp=M*ones(size(x,1),1);
% pdfexp(insidept==1) = evaluate_polyND(maxentpoly,x(insidept==1,:));
% pdfexp(insidept==0) = normpdfexp(x(insidept==0,:),m,P);

pdfexp=evaluate_polyND(maxentpoly,x);