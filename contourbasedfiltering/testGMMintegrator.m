close all
GMMfitter = GMMFitDataSet(dsX.X,dsX.p);
% GMM = GMMfitter.FitGMM_kmean_equalwt(5);
GMM = GMMfitter.FitGMM_kmeans_optimwt(5)


figure
GMMfitter.plotGMMpoints([1,2],'ro')
Xmc1 = random(MyGmm2MatlabGMM(GMM),1000);
hold on
plot(Xmc1(:,1),Xmc1(:,2),'b.')

for i=1:100
    Xmc1 = random(MyGmm2MatlabGMM(GMM),100*i);
    mean((pdfnorm.func(Xmc1))./GaussSumMix(Xmc1,GMM))
end

GMMhull = GMMfitter.SetGMM_kmeans_Hull();
figure
GMMfitter.plotGMMpointsHUll([1,2],3,'ro')
% Xmc1 = random(MyGmm2MatlabGMM(GMM),10000);
% hold on
% plot(Xmc1(:,1),Xmc1(:,2),'b.')

Z=[];
for i=1:GMM.Ngcomp
    Z=[Z;diag(sqrtm(GMM.Px{i}))'];
    cond(GMM.Px{i})
end


s = ComputeIntegral_importantGMM(GMM,pdfnorm.func,@(m,P)GH_pts(m,P,4))

normalize_exp_pdf(pdfnorm,dsX.X,dsX.p,mquad,Pquad,'GMM_MC')
