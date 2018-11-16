function f=constantfunc(x,const)
[r,c] = size(x);

if r==1 || c==1
    x=x(:)';
end

[r,c] = size(x);
f= const * ones(r,1);