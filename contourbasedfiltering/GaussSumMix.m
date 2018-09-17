function f = GaussSumMix(x,Ngcomp,wg,Mcomps,Pcomps)
[r,c]=size(x);
if r==1 || c==1
    x=x(:)';
else
    
end
[Np,dim]=size(x);

f=zeros(Np,1);
for i=1:Ngcomp
   f=f+wg(i)*mvnpdf(x,Mcomps{i}',Pcomps{i});
end