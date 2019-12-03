function h=radmodel_coe(xcoe,MU)

x = elm2rv_a(xcoe,0,MU);
x=x(1:3);
x=x(:);
% [xenu,~,~]=vec_radar_coordchange(x,PolarPosition,'ecef2local');
% x=x'-ecef_ref(Srad,:);
% xenu = ecef2enu(x*1e3,ecef_ref(Srad,:)*1e3)/1000;

% xenu=x-[6731,0,0]';
xenu=x;

xenu(:)';
r=norm(xenu);
th=atan2(xenu(1),xenu(2));
phi=atan2(sqrt(xenu(1)^2+xenu(2)^2),xenu(3));


h=[th;phi];

end