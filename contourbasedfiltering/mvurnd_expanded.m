function X = mvurnd_expanded(LB,UB,N,fac)
LB=LB(:);
UB=UB(:);
for i=1:length(LB)
   if LB(i)>0
       LB(i)=LB(i)*(1-fac);
   end
   if LB(i)<0
       LB(i)=LB(i)*(1+fac);
   end
   if UB(i)>0
       UB(i)=UB(i)*(1+fac);
   end
   if UB(i)<0
       UB(i)=UB(i)*(1-fac);
   end
end
% LB
% UB
dim = length(LB);
X = rand(N,dim);
X = X.*repmat((UB'-LB'),N,1)+repmat(LB',N,1);

end