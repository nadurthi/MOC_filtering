function [Y,yprobs]=filterX_inBox(X,probs,lb,ub)
[N,dim] = size(X);

lb = lb(:)';
ub = ub(:)';

y1 = sum(X>repmat(lb,N,1),2)==dim;
y2 = sum(X<repmat(ub,N,1),2)==dim;

Y= X(y1&y2,:);
if isempty(probs)
    yprobs = probs;
else
    yprobs = probs(y1&y2);
end

