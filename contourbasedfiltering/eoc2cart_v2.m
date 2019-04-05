function [rvec,vvec]=eoc2cart_v2(p,f,g,L,h,k,MU)

alpha2=h^2-k^2;
s2=1+h^2+k^2;
w=1+f*cos(L)+g*sin(L);

r=p/w;

rvec=[
(r/s2)*(cos(L)+alpha2*cos(L)+2*h*k*sin(L) )    

]