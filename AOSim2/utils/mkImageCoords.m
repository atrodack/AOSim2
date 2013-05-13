function [X1,X2,R] = mkImageCoords(I,SCALE,CENTER)

% function [X1,X2,R] = mkImageCoords(I,[SCALE],[CENTER])
%
% This function takes an image and an optional SCALE and CENTER (1,2)
% and returns a matching X, Y, and R coordinate overlay.
%
% 20061127: JLC: Finally wrote this stupid thing.
% 20081204: JLC: Changed definition of SCALE to be size of pixel.

x1 = 1:size(I,1);
x2 = 1:size(I,2);

if(nargin>2)
    x1 = x1 - CENTER(1);
    x2 = x2 - CENTER(2);
end

if(nargin>1)
%     x1 = x1/SCALE;
%     x2 = x2/SCALE;
    x1 = x1*SCALE;
    x2 = x2*SCALE;
end

[X2,X1] = meshgrid(x2,x1);


if(nargout==3)
    R = sqrt(X1.^2 + X2.^2);
end

return;
