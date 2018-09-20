function s = ComputeIntegral_importantGMM(GMM,func,genPointsSingleGauss)
% GMM.w;
% GMM.mx;
% GMM.Px;
% GMM.Ngcomp;

s=0;
for Nai = 1:GMM.Ngcomp
    alphai = GMM.w(Nai);
    [X,w]=genPointsSingleGauss(GMM.mx{Nai},GMM.Px{Nai});
    w=w(:);
    num1=func(X);
    dem1=GaussSumMix(X,GMM);
    s = s + alphai*sum(w.*(num1(:)./dem1(:)));
    
end

end