function [Y,probsY]=MeasUpdt_character(X,probs,z,model)
% the filter is implemented always using discrete - discrete models

probsXp = zeros(size(probs));
for i=1:size(X,1)
    probsXp(i) = model.z_pdf(z,X(i,:)')*probs(i);
end

%% Estimate normalizing constant

probsXp=probsXp/sum(probsXp);

[pdf,~] = get_interp_pdf(X,probsXp,4);

try
normpdf = normalize_exp_pdf(pdf,X,'dummyMC');
catch
    keyboard
end


%% Re-sample/ regenerate points

Y=X;
probsY=normpdf.func(Y);

% probsX=normpdf.func(X);
% [m,Px] = MeanCov(X,probsX/sum(probsX));
% Stdx=10*sqrtm(Px);
% [x1,x2]=meshgrid(linspace(m(1)-Stdx(1),m(1)+Stdx(1),model.Ngrid),linspace(m(2)-Stdx(2),m(2)+Stdx(2),model.Ngrid));
% 
% a=model.a;
% b=model.b;
% Y=[reshape(x1, a*b,1 ),reshape(x2, a*b,1 )];
% probsY=normpdf.func(Y);
