function X=gen_uniform_grid(n,dim)
wu=ones(n,1);
x=linspace(-1,1,n)';
X=[];
W=[];
for i=1:dim
[X,W]=tens_prod_vec(X,x,W,wu);
end