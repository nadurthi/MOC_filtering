function h=radmodel_2D(xnorm,constants)
xnorm=xnorm(:);

x = [xnorm(1:2)*constants.normX2trueX;xnorm(3:4)*constants.normV2trueV];
x=x(1:2);
x=x(:);
% [xenu,~,~]=vec_radar_coordchange(x,PolarPosition,'ecef2local');
% x=x'-ecef_ref(Srad,:);
% xenu = ecef2enu(x*1e3,ecef_ref(Srad,:)*1e3)/1000;

xenu=x;

% xenu(:)';
r=norm(xenu);
th=atan2(xenu(2),xenu(1));
% phi=atan2(sqrt(xenu(1)^2+xenu(2)^2),xenu(3));


h=[r;th];
% h=th;

end