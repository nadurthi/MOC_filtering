function f = GaussSumMix(x,GMM)
[r,c]=size(x);
if r==1 || c==1
    x=x(:)';
else
    
end
[Np,dim]=size(x);

f=zeros(Np,1);
for i=1:GMM.Ngcomp
    mm = GMM.mx{i};
    mm=mm(:);
   f=f+GMM.w(i)*mvnpdf(x,mm',GMM.Px{i});
end