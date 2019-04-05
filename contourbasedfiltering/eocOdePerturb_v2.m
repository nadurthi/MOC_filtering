function dx = eocOdePerturb_v2(x,MU,Re,J2)
% x=[p,f,g,L,h,k]
p=x(1);
f=x(2);
g=x(3);
L=x(4);
h=x(5);
k=x(6);

alpha2=h^2-k^2;
s2=1+h^2+k^2;
w=1+f*cos(L)+g*sin(L);

r=p/w;
sinphi = 2*(h*sin(L)-k*cos(L))/s2;

P2sinphi = 3/(1-sinphi^2)*(0.5*sinphi*(3*sinphi^2-1)-0.5*(5*sinphi^3-3*sinphi));

deltar=-3*MU*J2*Re^2/(2*r^4)*(1-12*(h*sin(L)-k*cos(L))^2/(1+h^2+k^2));
deltat=-12*MU*J2*Re^2/(r^4)*((h*sin(L)-k*cos(L))*(h*cos(L)+k*sin(L))/(1+h^2+k^2));
deltan=-6*MU*J2*Re^2/(r^4)*( ((1-h^2-k^2)*(h*sin(L)-k*cos(L)) )/(1+h^2+k^2));





dpdt = 2*p/w*sqrt(p/MU)*deltat;
dfdt = sqrt(p/MU)*(deltar*sin(L)+((w+1)*cos(L)+f)*deltat/w-(h*sin(L)-k*cos(L))*g*deltan/w   );
dgdt = sqrt(p/MU)*(-deltar*cos(L)+((w+1)*sin(L)+g)*deltat/w+(h*sin(L)-k*cos(L))*g*deltan/w   );

dhdt = sqrt(p/MU)*s2*deltan*cos(L)/(2*w);
dkdt = sqrt(p/MU)*s2*deltan*sin(L)/(2*w);
dLdt = sqrt(MU*p)*(w/p)^2+(1/w)*sqrt(p/MU)*(h*sin(L)-k*cos(L))*deltan;




dx=[dpdt,dfdt,dgdt,dLdt,dhdt,dkdt]';



