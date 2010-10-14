function I = airyPatternSecondary(theta,obs)

% function I = airyPatternSecondary(theta,obs) : The Airy pattern intensity.
%
% INPUT:
%   theta: the set of angels to evaluate measured in lambda/D.
%
% OUTPUT:
%   I: The intensity at theta normalized to the peak.
%
% Johanan L. Codona: 20051003:  Steward Observatory
% 20080926: obs!

theta = pi*(theta+eps);

Ap = besselj(1,theta)./(theta);
As = obs^2*besselj(1,theta*obs)./(theta*obs);

I = 4*(Ap-As).^2;
