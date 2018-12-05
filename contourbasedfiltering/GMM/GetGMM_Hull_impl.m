function GMMhull = GetGMM_Hull_impl(X,NgcompMax,method,mineigvalsetsqr)
[N,dim] = size(X);
switch nargin
    case 4
        mineigvalsqr=mineigvalsetsqr;
    otherwise
        mineigvalsqr=0.1;
end


% idx = GenerateClusterIndexes(X,NgcompMax,method);
idx = GenerateClusterIndexes_merging(X,NgcompMax,method);

Ngcomp = max(idx);
MF =cell(Ngcomp,1);
PF =cell(Ngcomp,1);
AF=cell(Ngcomp,1);
BF=cell(Ngcomp,1);

for i=1:Ngcomp
    xx=X(idx==i,:) ;
    Nc = size(xx,1);
    ww = ones(Nc,1)/Nc;
    [mcp,pcp]=MeanCov(xx,ww/sum(ww));
    xx = xx';
    [n,m] = size(xx);
    %     cvx_begin
    %         variable A(n,n) symmetric
    %         variable b(n)
    %         maximize( det_rootn( A ) )
    %         subject to
    %         norms( A * xx + b * ones( 1, m ), 2 ) <= 1;
    %     cvx_end
    [A,b]=computecvxCOVparallel(xx,n,m);
    
    % keyboard
    
    
    MF{i} = -inv(A)*b;
    
%     [u,s,v]=svd(A);
%     sd = diag(s);
%     sd(:)=5;
%     s=diag(sd);
%     A = u*s*v';
    
    PF{i} = inv(A*A);
        [u,s]=eig( PF{i} );
    sd = diag(s);
    sd(sd<mineigvalsqr)=mineigvalsqr;
    s=diag(sd);
    PF{i} = u*s*u';
    A =sqrtm(inv(PF{i}));
    b=-A*MF{i};
%     keyboard
    AF{i} = A;
    BF{i} = b;
end
GMMhull.w = ones(Ngcomp,1)/Ngcomp;
GMMhull.mx = MF;
GMMhull.Px = PF;
GMMhull.Ngcomp = Ngcomp;
GMMhull.A = AF;
GMMhull.b = BF;
GMMhull.idx=idx;
GMMhull.mineigvalsqr=mineigvalsqr;

% GMMhull=mergeGMMhulls(GMMhull,X);
GMMhull=addmoreGMMhulls(GMMhull,X);

end


function [A,b]=computecvxCOVparallel(xx,n,m)
cvx_begin
    variable A(n,n) symmetric
    variable b(n)
    maximize( det_rootn( A ) )
    subject to
    norms( A * xx + b * ones( 1, m ), 2 ) <= 1;
cvx_end

end