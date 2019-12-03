function xk1=processmodel_2body_2Dbackward(dt,MU,tk,xk)
% xk is at tk
% prop from tk to tk+dt
xk=xk(:);
[ r, v, Ehat ] = FnG(tk, tk-dt, [xk(1:2);0], [xk(3:4);0], MU);


xk1=[r(1:2);v(1:2)];

end