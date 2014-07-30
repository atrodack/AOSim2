function C = centroid(M)

% function C = centroid(M): First moment of M in 2-D.

M = squeeze(M);

C = [];

% M = M - mean(mean(M));
M1 = sum(M(:));
% M2 = mean(mean(abs(M).^2));

x1 = 1:size(M,1);
x2 = 1:size(M,2);

[X1,X2] = meshgrid(x1,x2);

C = [sum(sum(X1.*M)),sum(sum(X2.*M))]/M1;
