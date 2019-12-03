function f=treepredict_withbnds(X,tree,LB,UB)
LB=LB(:)';
UB=UB(:)';
[N,dim]=size(X);

y1 = sum(X>repmat(LB,N,1),2)==dim;
y2 = sum(X<repmat(UB,N,1),2)==dim;

indin = y1&y2;
indout = ~(y1&y2);

II=1:N;

Indin = II(indin);
Indout = II(indout);

Xin = X(indin,:);
Xout = X(indout,:);

fin = tree.predict(Xin);
fout = zeros(length(Indout),1);

f=zeros(N,1);
f(Indin)=fin;
f(Indout)=fout;

end