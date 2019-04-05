function [a,e,i,Om,om,M]=eoc2coe(p,f,g,L,h,k)

a=p/(1-f^2-g^2);
e=sqrt(f^2+g^2);
i=2*atan(sqrt(h^2+k^2));
Om=atan(k/h);
om=atan(g/f)-atan(k/h);
M=L-atan(g/f);









