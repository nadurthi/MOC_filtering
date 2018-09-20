function f=eval_truepdf_from_normpdf(xtrue,pdfnorm)
if any(size(xtrue)==1)
   xtrue=xtrue(:)'; 
end
[n,dim]=size(xtrue);
f=zeros(n,1);
for i=1:n
    xnorm=pdfnorm.trueX2normX(xtrue(i,:));
    pn = pdfnorm.func(xnorm);
    f(i) = pdfnorm.normprob2trueprob(pn);
end
