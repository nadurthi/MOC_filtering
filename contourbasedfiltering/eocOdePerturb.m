function dx = eocOdePerturb(x,MU,Re,J2)
p=x(1);
f=x(2);
g=x(3);
L=x(4);
h=x(5);
k=x(6);


s2=1+h^2+k^2;
w=1+f*cos(L)+g*sin(L);

r=p/w;
sinphi = 2*(h*sin(L)-k*cos(L))/s2;

P2sinphi = 3/(1-sinphi^2)*(0.5*sinphi*(3*sinphi^2-1)-0.5*(5*sinphi^3-3*sinphi));

dRdp = 3*MU/(2*w*r^2)*J2*(Re/r)^2*(12*(h*sin(L)-k*cos(L))^2/s2^2-1);
dRdf = -3*MU*cos(L)/(2*w*r)*J2*(Re/r)^2*(12*(h*sin(L)-k*cos(L))^2/s2^2-1);
dRdg = -3*MU*sin(L)/(2*w*r)*J2*(Re/r)^2*(12*(h*sin(L)-k*cos(L))^2/s2^2-1);

dRdL = -2*MU/(r*s2)*(h*cos(L)+k*sin(L))*J2*(Re/r)^2*P2sinphi-3*MU/(r*w)*(g*cos(L)-f*sin(L))*J2*(Re/r)^2*P2sinphi;

dRdh = -2*MU/(r*s2^2)*((1-h^2+k^2)*sin(L)+2*h*k*cos(L))*J2*(Re/r)^2*P2sinphi;
dRdk = 2*MU/(r*s2^2)*((1+h^2-k^2)*cos(L)+2*h*k*sin(L))*J2*(Re/r)^2*P2sinphi;







dpdt = 2*sqrt(p/MU)*(-g*dRdf+f*dRdg+dRdL);
dfdt = 1/sqrt(MU*p)*( 2*p*g*dRdp-(1-f^2-g^2)*dRdg-g*s2/2*(h*dRdh+kdRdk)+(f+(1+w)*cos(L))*dRdL  );
dgdt = 1/sqrt(MU*p)*( -2*p*f*dRdp+(1-f^2-g^2)*dRdf+f*s2/2*(h*dRdh+k*dRdk)+(g+(1+w)*sin(L))*dRdL );

dLdt = sqrt(MU*p)*(w/p)^2+s2/(2*sqrt(MU*p))*(h*dRdk+k*dRdk);
dhdt = s2/(2*sqrt(MU*p))*(h*(g*dRdf-f*dRdg-dRdL)-s2/2*dRdk  );
dkdt = s2/(2*sqrt(MU*p))*(k*(g*dRdf-f*dRdg-dRdL)+s2/2*dRdh  );



dx=[dpdt,dfdt,dgdt,dLdt,dhdt,dkdt]';



