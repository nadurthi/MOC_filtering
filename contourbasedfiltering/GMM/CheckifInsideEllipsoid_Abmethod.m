function ind = CheckifInsideEllipsoid_Abmethod(X,A,b,factor)
    [N,dim] = size(X);
    ind = zeros(N,1);
    for i=1:N
        ind(i)=norm(A*X(i,:)'+b,2)<=1*factor;
    end
end