function pest = estimateprobknn(x,X,p,Nn)
% keyboard
idx = knnsearch(X,x,'K',Nn);
if any(size(x)==1)
    pest = mean(p(idx));
else
    pest = mean(p(idx),2);
end
