function ind = CheckifInsideEllipsoid(X,m,P)
    m=m(:);
    invP  = inv(P);
    const = 1/sqrt(det(2*pi*P));
    [N,dim] = size(X);

    xt = sqrtm(P)*[1;zeros(dim-1,1)]+m;
    pt = const*exp(-0.5*(xt-m)'*invP*(xt-m));

    p = zeros(N,1);
    for i=1:N
        x = X(i,:)';
        p(i) = const*exp(-0.5*(x-m)'*invP*(x-m));
    end

    ind = zeros(N,1);
    ind(p>pt)=1;

end
