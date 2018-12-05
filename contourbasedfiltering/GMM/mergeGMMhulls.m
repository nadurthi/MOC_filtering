function GMMhull=mergeGMMhulls(GMMhull,X)
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

%first do single tons
flgalldone=0;
while(1)
    
    for i=1:GMMhull.Ngcomps
        Xmc=mvnrnd(GMMhull.mx{i},GMMhull.Px{i},Nmc);
        %check if these points are in any other of the hulls
        flginothergmm=0;
        for j=1:GMMhull.Ngcomps
            if i==j
                continue
            end
            A=GMMhull.A{i};
            b=GMMhull.b{i};
            
            S = sqrt( sum((A*Xmc'+repmat(b(:),1,Nmc)).^2,1) );
            if any(S<=1)
                flginothergmm=1;
                break
            end
        end
        if flginothergmm==1
            continue
        end
        % if not in other gmm then find the closest gmm
        % constru
    end
    
    if flgalldone==1
        break
    end
end

