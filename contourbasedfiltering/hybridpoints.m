function [X,w]=hybridpoints(mx,Px,pt1,fac1,pt2,fac2)

[X1,w1]=GH_points(mx,fac1^2*Px,pt1);
[X2,w2]=GH_points(mx,fac2^2*Px,pt2);

X=[X1;X2];
w=[w1;w2];