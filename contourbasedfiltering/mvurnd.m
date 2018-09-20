function X = mvurnd(LB,UB,N)
LB=LB(:);
UB=UB(:);

dim = length(LB);
X = rand(N,dim);
X = X.*repmat((UB'-LB'),N,1)+repmat(LB',N,1);

end