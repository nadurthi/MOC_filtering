function xk1=processmodel_2body_backward(dt,MU,tk,xk)
% xk is at tk
% prop from tk to tk+dt
xk=xk(:);
[ r, v, Ehat ] = FnG(tk, tk-dt, xk(1:3), xk(4:6), MU);


xk1=[r;v];

end