function OE = rv2elm_a(x,MU)

rbar=x(1:3);
vbar=x(4:6);

rbar=rbar(:);
vbar=vbar(:);

r=norm(rbar);
v=norm(vbar);

khat = [0,0,1]';

hbar = cross(rbar,vbar);
h = norm(hbar);

ebar = (1/MU)*( (v^2-MU/r)*rbar-dot(rbar,vbar)*vbar );
e=norm(ebar);

nbar = cross(khat,hbar);
n = norm(nbar);

p=h^2/MU;
a=p/(1-e^2);

inc = acos(dot(hbar,khat)/h);

Om=acos(nbar(1)/n);
if nbar(2)<0
    Om=2*pi-Om;
end

om=acos(dot(nbar,ebar)/(n*e));
if ebar(3)<0
  om=2*pi-om;
end

f = acos(dot(ebar,rbar)/(e*r));
if dot(rbar,vbar)<=0
    f=2*pi-f;
end

E = acos( (e+cos(f))/(1+e*cos(f)) );
if f>pi
    E=2*pi-E;
end

M=E-e*sin(E);



OE=[a,e,inc,Om,om,M];


