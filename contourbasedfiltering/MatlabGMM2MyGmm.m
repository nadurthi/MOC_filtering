function GMM = MatlabGMM2MyGmm(gmm)
GMM.w = gmm.ComponentProportion';
GMM.Ngcomp = gmm.NumComponents';

dim = gmm.NumVariables;

GMM.mx = cell(GMM.Ngcomp,1);
GMM.Px = cell(GMM.Ngcomp,1);
for i=1:GMM.Ngcomp
    GMM.mx{i} = gmm.mu(i,:);
    GMM.Px{i} = gmm.Sigma(:,:,i);
end
