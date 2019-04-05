function Xto = boxShift_working(Xfrom,fromlb,fromub,tolb,toub)
% trasnform box from box limits --> to box limits
[r,c]=size(Xfrom);
if r==1 || c==1
   Xfrom=Xfrom(:)'; 
end

maxX = fromub;
minX = fromlb;
maxX=maxX(:);
minX=minX(:);

LB=tolb(:);
UB=toub(:);
A = diag((UB-LB)./(maxX-minX));
m = -(minX.*(UB-LB))./(maxX-minX)+LB;

for i=1:size(Xfrom,1)
    Xfrom(i,:) = A*Xfrom(i,:)'+m(:);
end

Xto = Xfrom;
% scale to size of tobox

