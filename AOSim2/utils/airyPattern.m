function I = airyPattern(theta)

% function I = airyPattern(theta) : The Airy pattern intensity.
%
% INPUT:
%   theta: the set of angels to evaluate measured in lambda/D.
%
% OUTPUT:
%   I: The intensity at theta normalized to the peak.
%
% Johanan L. Codona: 20051003:  Steward Observatory

theta = pi*(theta+eps);

I = 4*(besselj(1,theta+eps)./(theta+eps)).^2;
