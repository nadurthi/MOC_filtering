function ind = CheckifInsideEllipsoid_Abmethod(X,A,b,factor)
    [N,dim] = size(X);
    ind = zeros(N,1);
%     for i=1:N
%         ind(i)=norm(A*X(i,:)'+b,2)<=1*factor;
%     end
% try
    S = sqrt( sum((A*X'+repmat(b(:),1,N)).^2,1) );
% catch
%     keyboard
% end

    ind(S<=1*factor)=1;
end