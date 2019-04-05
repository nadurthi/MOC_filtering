function h=radmodel_orb2angles(xk,MU)
% [a,e,i,om,Om,M]
a=xk(1);
e = xk(2);
inc = xk(3);
om = xk(4);
Om = xk(5);
Mprev = xk(6);
[r0,v0] = elm2rv(a,e,inc,Om,om,Mprev,0,MU);

% x=x(1:3);
x=r0(:);
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