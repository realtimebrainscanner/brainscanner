function [X] = center(X)
n = size(X,1);
mu = mean(X);
X = X - ones(n,1)*mu;
