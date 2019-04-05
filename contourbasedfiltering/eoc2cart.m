function [r,v]=eoc2cart(p,f,g,L,h,k,MU)

a=p/(1-f^2-g^2);
e=sqrt(f^2+g^2);
inc=2*atan(sqrt(h^2+k^2));
Om=atan(k/h);
om=atan(g/f)-atan(k/h);
M=L-atan(g/f);

[r,v] = elm2rv(a,e,inc,Om,om,M,0,MU);