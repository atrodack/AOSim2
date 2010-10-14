function V = demean(V)

% function V = demean(V)
%   remove the mean. 
%  JLC: 20070531

V = V - mean(V(:));

return;
