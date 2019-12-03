function xk1=processmodel_2bodycoe(dt,MU,tk,xk)
% xk is at tk
% prop from tk to tk+dt
% xk is coe
% xk=xk(:);
% ek = xk(2);
% Mk = xk(6);
xkcart = elm2rv_a(xk,0,MU);

[ r, v, Ehat ] = FnG(tk, tk+dt, xkcart(1:3), xkcart(4:6), MU);

xk1 = rv2elm_a([r(:)',v(:)'],MU);

end