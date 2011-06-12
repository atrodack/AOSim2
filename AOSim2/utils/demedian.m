function V = demedian(V)

% function V = demean(V)
%   remove the mean. 
%  JLC: 20070531

V = V - median(V(:));

return;
