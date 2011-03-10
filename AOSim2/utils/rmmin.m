function V = rmmin(V)

% function V = rmmin(V)
%   remove the min value so the smallest value is 0.0. 
%  JLC: 20110202

V = V - min(V(:));

return;
