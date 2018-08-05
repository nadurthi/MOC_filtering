clc
clear
close all
digits(16)
mxentpoly_x =[1     7   0     0;
     15     3     2     1;
     900     1     3     2;
     4     2    1     3];
%  A=2*eye(3);
%  A=diag([1;2;3]);

A=randn(3).^2;
A=A*A';
A=[1.76325784693713          3.99092800170445          2.90979133952959
          3.99092800170445          9.14154078315538          7.32655941718451
          2.90979133952959          7.32655941718451          10.0599529650316];

% A=100*A;
% A=round(A);

A=diag(diag(A));
A(1,2)=40;
A(2,1)=40;
% A(1,3)=120.340;
% A(3,1)=-120.340;
A(2,3)=-0.00120340;
A(3,2)=-0.00120340;
A
eig(A)
mu=[1;3;5];

mxentpoly_y=linear_transform_poly(mxentpoly_x,inv(A),-inv(A)*mu(:));
mxentpoly_x2=linear_transform_poly(mxentpoly_y,A,mu(:));

disp('---++++---')
Xtest1=[2.3;3;7];
evaluate_polyND(mxentpoly_x,Xtest1)
evaluate_polyND(mxentpoly_y,A*Xtest1+mu(:))

disp('-------')
Xtest1=[5;10.5;14.7];
evaluate_polyND(mxentpoly_x,Xtest1)
evaluate_polyND(mxentpoly_x2,Xtest1)



%% 
% 
% Ainv = inv(A);
% a11=Ainv(1,1);
% a12=Ainv(1,2);
% a21=Ainv(2,1);
% a22=Ainv(2,2);
% a33=Ainv(3,3);
% 
% 
% Py=[nchoosek(7,0)*a11^7*a12^0,      7,0,0;
%     nchoosek(7,1)*a11^6*a12^1,      6,1,0;
%     nchoosek(7,2)*a11^5*a12^2,      5,2,0;
%     nchoosek(7,3)*a11^4*a12^3,      4,3,0;
%     nchoosek(7,4)*a11^3*a12^4,      3,4,0;
%     nchoosek(7,5)*a11^2*a12^5,      2,5,0;
%     nchoosek(7,6)*a11^1*a12^6,      1,6,0;
%     nchoosek(7,7)*a11^0*a12^7,      0,7,0;
% ];
% 
% Xtest1=[1;1;1];
% evaluate_polyND(mxentpoly_x,Xtest1)
% evaluate_polyND(mxentpoly_y,A*Xtest1+mu(:))
% evaluate_polyND(Py,A*Xtest1+mu(:))
% 
% mxentpoly_x2=linear_transform_poly(Py,A,mu(:));
