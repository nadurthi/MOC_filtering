function f = evalKnnMean_bykdtree(Mdl,Ytrain,Xtest)
[r,c]=size(Xtest);

if r==1 || c==1
    Xtest = Xtest(:)';
end
[Ntest,dimtest]=size(Xtest);
Nn=dimtest+1;


idx = knnsearch(Mdl,Xtest,'K',Nn);

Xtrain=Mdl.X;

f=zeros(Ntest,1);

for i=1:Ntest
    dists = sqrt(sum((Xtrain(idx(i,:),:)-repmat(Xtest(i,:),Nn,1)).^2,2));
    wts = 1./dists;
    wts = wts/sum(wts);
    f(i) = sum(wts.*Ytrain(idx(i,:)));
end