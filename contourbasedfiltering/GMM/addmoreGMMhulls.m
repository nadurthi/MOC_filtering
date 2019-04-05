function GMMhull=addmoreGMMhulls(GMMhull,X)
mineigvalsqr = GMMhull.mineigvalsqr;

% join all the GMM that do not overlap
% GMMhull.w = ones(Ngcomp,1)/Ngcomp;
% GMMhull.mx = MF;
% GMMhull.Px = PF;
% GMMhull.Ngcomp = Ngcomp;
% GMMhull.A = AF;
% GMMhull.b = BF;

[N,dim]=size(X);
if dim==2
    Nmc=500;
end
if dim==4
    Nmc=1000;
end
if dim==6
    Nmc=10000;
end

Ng= GMMhull.Ngcomp;

MF={};
PF={};
AF={};
BF={};
k=1;
D=zeros(GMMhull.Ngcomp,GMMhull.Ngcomp);

for i=1:GMMhull.Ngcomp
%     Xmc=mvnrnd(GMMhull.mx{i},GMMhull.Px{i},Nmc);
%     %check if these points are in any other of the hulls
%     flginothergmm=0;
    for j=i+1:GMMhull.Ngcomp
        [i,j,GMMhull.Ngcomp]
        A1=GMMhull.A{i};
        A2=GMMhull.A{j};
        b1=GMMhull.b{i};
        b2=GMMhull.b{j};

        cvx_begin
            variables x1(dim) x2(dim)
            minimize( norm(x1-x2) )
            subject to
            norm( A1 * x1 + b1, 2 ) <= 1;
            norm( A2 * x2 + b2, 2 ) <= 1;
        cvx_end
        D(i,j)=norm(x1-x2);
        D(j,i)=D(i,j);
    end
end
% now get the clusters using D
Dd=D;
Dd(Dd<=0.01)=1;
Dd((Dd>0.01) & (Dd~=1))=0;


ClusterBags={};
k=1;
II=1:GMMhull.Ngcomp;
for i=1:GMMhull.Ngcomp
   ClusterBags{k}=II(Dd(i,:)==1);
   k=k+1;
end


while(1)
    Delll=[];
    for i=1:1:length(ClusterBags)
        for j=i+1:length(ClusterBags)
            C = intersect(ClusterBags{i},ClusterBags{j});
            if isempty(C)==0
                ClusterBags{i} = unique([ClusterBags{i},ClusterBags{j}]);
                ClusterBags{j}=[];
                Delll=horzcat(Delll,j);
            end
        end
    end
    if isempty(Delll)==1
        break
    else
        ClusterBags(Delll)=[];
    end
end

k=1;
while(1)
    if length(ClusterBags)==1
        break
    end
    flgbrk=0;
    for i=1:1:length(ClusterBags)
        for j=i+1:length(ClusterBags)
            gmmsi = ClusterBags{i};
            gmmsj = ClusterBags{j};
            mind=1000;
            mini=0;
            minj=0;
            for si=1:length(gmmsi)
                for sj=1:length(gmmsj)
                    mxi =GMMhull.mx{gmmsi(si)};
                    mxj =GMMhull.mx{gmmsj(sj)};
                    d=norm(mxi-mxj);
                    if d<mind
                        mind=d;
                        mini=gmmsi(si);
                        minj=gmmsj(sj);
                    end
                end
            end
%             keyboard
            A1=GMMhull.A{mini};
            A2=GMMhull.A{minj};
            b1=GMMhull.b{mini};
            b2=GMMhull.b{minj};

            S1 = sqrt( sum((A1*X'+repmat(b1(:),1,N)).^2,1) );
            S2 = sqrt( sum((A2*X'+repmat(b2(:),1,N)).^2,1) );
            Y=X((S1<=1.1) | (S2<=1.1),:);
            Y=Y';
            [n,m] = size(Y);

            [A,b]=computecvxCOVparallel(Y,n,m);
            MF{k,1} = -inv(A)*b;

            PF{k,1} = inv(A*A);
            [u,s]=eig( PF{k} );
            sd = diag(s);
            sd(sd<mineigvalsqr)=mineigvalsqr;
            s=diag(sd);
            PF{k,1} = u*s*u';
            A =sqrtm(inv(PF{k,1}));
            b=-A*MF{k};
            AF{k,1} = A;
            BF{k,1} = b;
            k=k+1;
            
            ClusterBags{i}=unique([ClusterBags{i},ClusterBags{j}]);
            ClusterBags(j)=[];
            flgbrk=1;
            break
        end
        if flgbrk==1
            break
        end
    end
    
end
% keyboard

GMMhull.mx = vertcat(GMMhull.mx,MF);
GMMhull.Px = vertcat(GMMhull.Px,PF);
GMMhull.A = vertcat(GMMhull.A,AF);
GMMhull.b = vertcat(GMMhull.b,BF);
GMMhull.Ngcomp = length(GMMhull.mx);
GMMhull.w = ones(GMMhull.Ngcomp,1)/GMMhull.Ngcomp;

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
