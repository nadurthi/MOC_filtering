function [m,P]=GMMcell2array(GMM)
dim = length(GMM.mx{1});
m=zeros(GMM.Ngcomp,dim);
P=zeros(GMM.Ngcomp,dim^2);
for i=1:GMM.Ngcomp
    m(i,:)=GMM.mx{i};
    P(i,:)=reshape(GMM.Px{i},1,dim^2);
end