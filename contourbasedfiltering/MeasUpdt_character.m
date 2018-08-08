function [Y,probsY]=MeasUpdt_character(X,probs,mquad,Pquad,Nm,Tk,z,model)
% the filter is implemented always using discrete - discrete models

probsXp = zeros(size(probs));
for i=1:size(X,1)
    probsXp(i) = model.z_pdf(z,X(i,:)')*probs(i);
end

%% Estimate normalizing constant

[pdfnorm,pdftransF] = get_interp_pdf_0I(X,probs,mquad,Pquad,Nm,Tk,[])
y=pdftransF.trueX2normY(X);
py=pdfnorm.func(y);
probsY=pdftransF.normprob2trueprob(py);
%% Re-sample/ regenerate points

[Y,w] = GH_points(mquad,0.5^2*Pquad,5);

probsY=normpdf.func(Y);


