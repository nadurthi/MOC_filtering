function X1 = sphere4Dm(n)

N=n^3;
X1 = zeros(N,4);

% ph1,ph2,ph3,ph4,ph5
ph1 =linspace(0*pi/180,180*pi/180,n);
ph2 =linspace(0*pi/180,180*pi/180,n);
ph3 =linspace(0*pi/180,350*pi/180,n);

a=1;
b=1;


for k=1:length(ph1)
    for p=1:length(ph2)
        for q=1:length(ph3)
            X1(a,1)=cos(ph1(k));
            X1(a,2)=sin(ph1(k))*cos(ph2(p));
            X1(a,3)=sin(ph1(k))*sin(ph2(p))*cos(ph3(q));
            X1(a,4)=sin(ph1(k))*sin(ph2(p))*sin(ph3(q));
            a=a+1;
        end
    end
    
    %     keyboard
    
    
end


% H=10^6*X1(:,1)+10^5*X1(:,2)+10^4*X1(:,3)+10^3*X1(:,5)+10^2*X1(:,5)+X1(:,6);
% Hunq=unique(H);

X1=unique(X1,'rows');








end