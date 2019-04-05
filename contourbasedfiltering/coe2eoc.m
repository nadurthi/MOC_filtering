function [p,f,g,L,h,k]=coe2eoc(a,e,i,Om,om,M)

p=a*(1-e^2);
f=e*cos(Om+om);
g=e*sin(Om+om);
L=M+Om+om;
h=tan(i/2)*cos(Om);
k=tan(i/2)*sin(Om);






