function [Y,probs]=propagate_character_back(X,probs,dt,tk,model)
% propogate from tk to tk-dt using model.fback

Y=zeros(size(X));
probsy = zeros(size(probs));
for i=1:size(X,1)
    Y(i,:) = model.fback(dt,tk,X(i,:)');
    probsy(i) = probs(i);
end
