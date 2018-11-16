function xk1=processmodel_2body_orb(dt,MU,tk,xk)
% xk is at tk [a,e,i,om,Om,M]
% prop from tk to tk+dt
a=xk(1);
e = xk(2);
inc = xk(3);
om = xk(4);
Om = xk(5);
Mprev = xk(6);
[r0,v0] = elm2rv(a,e,inc,Om,om,Mprev,0,MU);
[ r, v, Ehat ] = FnG(tk, tk+dt, r0, v0, MU);

OE = cart2orbelem([r;v],MU);

xk1 = [OE.a,OE.e,OE.i,OE.om,OE.Om,OE.M];

end