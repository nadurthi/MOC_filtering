digits(50);
dsXX = DataSet(vpa(Xrepk1(1:50,:)),vpa(prepk1(1:50,:)),'TrueState');

% [mm,PP]=MeanCov(vpa(Xrepk1),vpa(prepk1/sum(prepk1)));
dsXX.AddMeanCov_to_OI_Trasform(vpa(mXk1),vpa(4^2*PXk1));

[Y,~]=dsXX.ApplyAffineTransform_Final2Original(dsXX.X(1:50,:),ones(50,1));

Xrepk1(1:50,:)-Y
%%
dsXX = DataSet(Xrepk1,prepk1,'TrueState');

% [mm,PP]=MeanCov(vpa(Xrepk1),vpa(prepk1/sum(prepk1)));
dsXX.AddMeanCov_to_OI_Trasform(mXk1,4^2*PXk1);

Y=dsXX.ApplyAffineTransform_Final2Original(dsXX.X,ones(50,1));

Xrepk1(1:50,:)-Y(1:50,:)
%%
XX=Xrepk1(1:50,:);
A=sqrtm(PP);
mm=mm(:)';
XX0I = zeros(size(XX));
for i=1:50
    XX0I(i,:) = A\(XX(i,:)-mm)';
end