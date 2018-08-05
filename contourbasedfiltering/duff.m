function dx=duff(t,x)
dx=zeros(2,1);
dx(1)=x(2);
dx(2)=-x(1)-1*x(1)^3;

end