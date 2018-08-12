Xmc=XMC(:,:,33);
figure
plot(Xmc(:,1),Xmc(:,2),'bo')

[mx,Px] = MeanCov(Xmc,ones(Nmc,1)/Nmc);

Xnmc = zeros(size(Xmc));
Psqrt_inv = sqrtm(inv(Px)); 
for i=1:Nmc
    Xnmc(i,:) = Psqrt_inv*(Xmc(i,:)'-mx(:));
end
figure
plot(Xnmc(:,1),Xnmc(:,2),'bo')

Xrthmc = zeros(size(Xmc));
for i=1:Nmc
    Xrthmc(i,1) = sqrt(sum(Xmc(i,:).^2));
    Xrthmc(i,2) = atan(Xmc(i,2)/Xmc(i,1));
end
figure
plot(Xrthmc(:,1),Xrthmc(:,2),'bo')


