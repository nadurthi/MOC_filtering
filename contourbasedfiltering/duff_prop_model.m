function x1=duff_prop_model(dt,x)
[t,x]=ode45(@duff,[0,dt],x);

x1=x(end,:)';