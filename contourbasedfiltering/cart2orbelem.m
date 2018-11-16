function OE = cart2orbelem(XX,MU)
XX=XX(:);
X=XX(1:3);
V=XX(4:6);

%% Input - X,V position and velocity of the point mass... 

rm = sqrt(X'*X); vm = sqrt(V'*V);

a = 1/(2/rm - vm^2/MU);

hbar = cross(X,V);
cbar = cross(V,hbar)-MU*X/rm;

e = sqrt(cbar'*cbar)/MU;

hm = sqrt(hbar'*hbar); ih = hbar/hm;

ie = cbar/(MU*e);

ip = cross(ih,ie);

i = acos(ih(3,1));

w = atan2(ie(3,1),ip(3,1));

Om = atan2(ih(1,1),-ih(2,1));

f=acos(dot(ip,X)/(1*rm));
if dot(V,X)<=0
    f=2*pi-f;
end

sig = X'*V/sqrt(MU);
% keyboard
E = atan2(sig/sqrt(a),1-rm/a);

M = E-e*sin(E);

OE.a=a;
OE.e=e;
OE.om=w;
OE.i=i;
OE.Om=Om;
OE.f=f;
OE.M=M;
OE.E=E;

% OE=[a,e,E,M,w,i,Om];

