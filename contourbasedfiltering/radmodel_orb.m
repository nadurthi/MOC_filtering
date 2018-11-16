function h=radmodel_orb(xorb,MU)
% xk is at tk [a,e,i,om,Om,M]
a=xorb(1);
e = xorb(2);
inc = xorb(3);
om = xorb(4);
Om = xorb(5);
M = xorb(6);
[r,v] = elm2rv(a,e,inc,Om,om,M,0,MU);


x=r;
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