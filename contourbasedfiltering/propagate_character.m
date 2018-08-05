function [Y,probs]=propagate_character(X,probs,dt,tk,model)
% the filter is implemented always using discrete - discrete models

Y=zeros(size(X));
probsy = zeros(size(probs));
for i=1:size(X,1)
    Y(i,:) = model.f(dt,tk,X(i,:)');
    probsy(i) = probs(i);
end
