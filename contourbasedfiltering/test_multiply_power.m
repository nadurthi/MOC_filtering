clc
clear
mxentpoly_x =[
     1     3     0     0
    15     3     2     1
   900     1     3     2
     4     2     1     3
     ];
 
 pw=3;
 testpt=[10,12,4];
 tic
 p1=tothepowerof_polyND(mxentpoly_x,pw) ;
S=binomexp_NDpoly(mxentpoly_x(2:end,:),mxentpoly_x(1,:),pw);
x=evaluate_polyND(mxentpoly_x,testpt);
x^pw
x^pw-vpa(evaluate_polyND(p1,sym(testpt)))
x^pw-evaluate_polyND2(p1,testpt)
x^pw-evaluate_polyND(S,testpt)
toc

%%
% clc
% v=sym('v',[1,3]);
% f=vpa(evaluate_polyND_sym(sym(p1),v));
% vpa(x^pw)
% subs(vpa(f),v,sym(testpt))
% vpa(evaluate_polyND(S,testpt))