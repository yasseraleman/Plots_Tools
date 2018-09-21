function [b,c] = my_linearfit(X,Y);
X = X(:);Y = Y(:);
n= size(X,1);
epxy = dot(X,Y);
epxx = dot(X,X);
sxx =sum(X)^2;
b = (n*epxy-sum(X)*sum(Y))/(n*epxx-sxx);
c = (epxx*sum(Y)-sum(X)*epxy)/(n*epxx-sxx); 
return;
