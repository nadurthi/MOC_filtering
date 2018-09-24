function [X,w]=mvnrnd_modf(m,P,N)
X = mvnrnd(m(:)',P,N);
w=ones(N,1)/N;


end