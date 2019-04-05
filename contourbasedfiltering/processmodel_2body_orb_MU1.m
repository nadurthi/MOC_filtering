function xk1=processmodel_2body_orb_MU1(dt,MU,RE,DU,VU,TU,tk,xk)
% xk is at tk [a,e,i,om,Om,M], a is in canonical units
% prop from tk to tk+dt
a=xk(1);
e = xk(2);
inc = xk(3);
om = xk(4);
Om = xk(5);
Mprev = xk(6);
[r0,v0] = elm2rv(a,e,inc,Om,om,Mprev,0,MU);
% [ r, v, Ehat ] = FnG(tk, tk+dt, r0, v0, MU);

opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,x]=ode45(@(t,x)twoBody_j2(t,x,MU,RE),[tk*TU,(tk+dt)*TU],[r0*DU,v0*VU],opts);

xend = x(end,:)';
xend(1:3)=xend(1:3)/DU;
xend(4:6)=xend(4:6)/VU;
OE = cart2orbelem(xend,MU);

xk1 = [OE.a,OE.e,OE.i,OE.om,OE.Om,OE.M];

end