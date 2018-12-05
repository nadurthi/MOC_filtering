function f = sensmodel2vec(x,hfunc,hn)
[r,c] = size(x);

if r==1 || c == 1
    x=x(:)';
end

[N,dim] = size(x);

f= zeros(N,hn);
for i=1:N
   f(i,:) = hfunc(x(i,:)'); 
end