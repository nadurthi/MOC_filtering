function y=fixedcoordpdfexp(normpdffunc,fixedstates,evalxs,fixedx,dim)
[r,c]=size(evalxs);
if r==1||c==1
    evalxs=evalxs(:)';
end

[N,~]=size(evalxs);

y=zeros(N,1);

x=zeros(N,dim);
x(:,fixedstates)=repmat(fixedx(:)',N,1);

evalstates = setdiff(1:dim,fixedstates);
x(:,evalstates) = evalxs;

y = normpdffunc(x);




