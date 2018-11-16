function GMM=getInitialGMM(x0,P0,N1d)
dim = length(x0);

sqP0 = sqrtm(P0);


if rem(N1d,2)==0
    N1d=N1d+1;
end


[X,w]=GH_points(0,1,N1d);
% X=2*X/max(abs(X));
% w=ones(N1d,1)/N1d;

for n0=1:length(w)
    if X(n0)==0
        break
    end
end

while(1)
    c2=0;
    for i=1:length(w)
        c2=c2+w(i)*X(i)^2;
    end
    c2=1-c2;
    
    if c2<0
        w(n0)=w(n0)+0.01;
        w=w/sum(w);
    else
        break
    end
    
end


F=[w(n0),2*w(n0+1:end)'];

A=[];
z=zeros(1,n0);
for i=1:1:n0-1
    d=z;
    d(i)=-1;
    d(i+1)=1;
    A=vertcat(A,d);
end
% x=fmincon(@(x)(3*w(n0)*x(1)^2+2*sum(w(n0+1:end).*X(n0+1:end).^4)+12*sum(w(n0+1:end).*(X(n0+1:end).^2).*x(2:3))+6*sum(w(n0+1:end).*(x(2:3).^2)) -4)^2,ones(n0,1),A,zeros(n0-1,1),F,c2,zeros(n0,1),200*ones(n0,1));

x=linspace(-5,5,100);
ytrue = normpdf(x);
s=fmincon(@(s)sum((gmm1d(x,w,X,s)-ytrue).^2),ones(n0,1),[],[],[],[],zeros(n0,1),200*ones(n0,1));




sigmas=zeros(length(w),1);
sigmas(n0:end)=s;
sigmas(1:n0-1)=s(end:-1:2);

w
wu=w;
u=X;
Xall=X;
Wall=w;
Sall=sigmas;
G=w;
for i=1:dim-1 
    [Xall,Wall]=tens_prod_vec(Xall,X,Wall,w);
    [Sall,G]=tens_prod_vec(Sall,sigmas,G,w);
    
end

GMM.w=Wall;
GMM.Ngcomp=length(Wall);
GMM.mx=cell(length(Wall),1);
GMM.Px=cell(length(Wall),1);
for i=1:length(Wall)
    GMM.mx{i}=Xall(i,:)'+x0;
    GMM.Px{i}=sqP0*diag(Sall(i,:))*sqP0';
end
sum(GMM.w)

% yest=gmm1d(x,w,X,s);
% figure
% plot(x,ytrue,'b',x,yest,'r');


% figure
% hold on
% states=[1,2];
% for i=1:GMM.Ngcomp
%     plot_nsigellip(GMM.mx{i}(states),GMM.Px{i}(states,states),1,'b',1)
% end
% axis equal
% axis square
% 
% Xt=mvnrnd(x0(:)',P0,10000);
% fest=evalGMM(GMM,Xt);
% ftrue = mvnpdf(Xt,x0(:)',P0);
% 
% figure
% plot3(Xt(:,1),Xt(:,2),ftrue,'bo')
% hold on
% plot3(Xt(:,1),Xt(:,2),fest,'r+')

end

function y=gmm1d(x,w,mm,sig2)
n0=length(sig2);
sigmas2  = zeros(1+2*(n0-1),1);
sigmas2(n0:end)=sig2;
sigmas2(1:n0-1)=sig2(end:-1:2);
y=0;
for i=1:length(w)
    y=y+w(i)*normpdf(x,mm(i),sigmas2(i));
end

end

