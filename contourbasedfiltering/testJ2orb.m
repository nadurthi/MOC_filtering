clc
close all
clear

format longg

digits(50)
%% constants
constants.radii=[6378.137,6378.137,6378.137];

constants.mu      = 3.986004418e5;     % Gravitational Const
constants.Re      = constants.radii(1);          % Earth radius (km)

constants.g0   = 9.8065;            % Sea-level acceleration, (m/s^2)

% Canonical Units
constants.muCan   = 1;
constants.RU      = constants.Re;
constants.TU      = sqrt(constants.RU^3 / constants.mu);
constants.VU      = constants.RU/constants.TU;

constants.trueA2normA=(constants.TU^2/constants.RU);
constants.normA2trueA=(constants.RU/constants.TU^2);

constants.trueV2normV=(constants.TU/constants.RU);
constants.normV2trueV=(constants.RU/constants.TU);

constants.trueX2normX=(1/constants.RU);
constants.normX2trueX=(constants.RU);

constants.trueT2normT=(1/constants.TU);
constants.normT2trueT=(constants.TU);

%%
x0true=[10000,10,3000,-0.1,7.4,2]';
OE0true = cart2orbelem(x0true,constants.mu)
tktrue = 24*3600;

%% prop true cart
RE=constants.Re;
MU=constants.mu;

opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,x]=ode45(@(t,x)twoBody_j2(t,x,MU,RE),[0,tktrue],x0true,opts);

%% prop true cart

opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,x]=ode45(@(t,x)twoBody_j2(t,x,MU,RE),[0,tktrue*constants.TU],x0true,opts);




[r0,v0] = elm2rv(a,e,inc,Om,om,Mprev,0,MU);
OE = cart2orbelem(xend,MU);


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
