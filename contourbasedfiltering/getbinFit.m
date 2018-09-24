function BinFit = getbinFit(X,probs,LB,UB,Nbin)
[N,dim]=size(X);

Xbins=gen_uniform_grid(Nbin,dim);
Xbins=Xbins+1;
Xbins=Xbins/2;
Xbins = Xbins.*repmat((UB'-LB'),N,1)+repmat(LB',N,1);

%% get avg estimate in each bin
for 
