figure
plot(XMC(:,1,23),XMC(:,2,23),'ro')
[m,P]=MeanCov(XMC(:,:,23),ones(Nmc,1)/Nmc);
cond(P)
Xn=zeros(Nmc,4);


Psqrt=sqrtm(P);
Psqrt_inv=inv(sqrtm(P));
for i=1:Nmc
    Xn(i,:) = Psqrt_inv*(XMC(i,:,23)-m(:)')';
    
end
figure
plot(Xn(:,1),Xn(:,2),'ro')

%% 
Xpmc = zeros(Nmc,2);
[m,P]=MeanCov(Xpmc,ones(Nmc,1)/Nmc);
m2 = Xpmc(knnsearch(Xpmc,m'),:);

P=zeros(2);
for i=1:Nmc
   P=P+1/Nmc * (Xpmc(i,:)-m2(:)')'*(Xpmc(i,:)-m2(:)'); 
end
figure
plot(Xpmc(:,1),Xpmc(:,2),'ro',m2(1),m2(2),'b+')
hold on
plot_nsigellip(m2,P,1,'b',2)


Xn=zeros(Nmc,2);
% Psqrt=sqrtm(P);
% PP = diag(diag(P));
PP = P;
Psqrt_inv=inv(sqrtm(1*PP));
for i=1:Nmc
    Xn(i,:) = Psqrt_inv*(Xpmc(i,:)-m2(:)')';
    
end


figure
plot(Xn(:,1),Xn(:,2),'ro')
hold on
plot_nsigellip(zeros(2,1),eye(2),1,'b',2)
hold off
axis equal
axis square
% axis([-4,4,-4,4])
% 