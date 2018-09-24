function gmm = MyGmm2MatlabGMM(GMM)
dim = length(GMM.mx{1});
P = GMM.w(:)';
Mu=zeros(GMM.Ngcomp,dim);
Sigma = zeros(dim,dim,GMM.Ngcomp);
for i=1:GMM.Ngcomp
    Mu(i,:) = GMM.mx{i};
    Sigma(:,:,i) = GMM.Px{i};
end
% Mu = [1 2;-3 -5];
% Sigma = cat(3,[2 0;0 .5],[1 0;0 1]);
% 
gmm = gmdistribution(Mu,Sigma,P);