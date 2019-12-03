function [X,pts1D,interpPoly1D]=getsparsePts_scaled(q,d,scale)

% d=2;
% dim=d;
% q=8;
pts1D=cell(1,q-d+2);
interpPoly1D = cell(1,q-d+2);

% digits(50)



for i=1:q-d+2
    if i==1
        pts1D{i}=scale*chebychevextrema(1);
    else
        pts1D{i}=scale*chebychevextrema(2^(i-1)+1);
    end
    interpPoly1D{i}=LagrangePoly1D_matpoly( pts1D{i} ,1,1);  % vpa(pts1D{i})  pts1D{i}
end
X=sparseGridConstructor(q,d,pts1D);
