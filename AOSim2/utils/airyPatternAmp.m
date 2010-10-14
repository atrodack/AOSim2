function I = airyPatternAmp(theta)

% function I = airyPatternAmp(theta) : The Airy pattern AMPLITUDE.
%
% INPUT:
%   theta: the set of angels to evaluate measured in lambda/D.
%
% OUTPUT:
%   I: The AMPLITUDE at theta normalized to the peak.
%
% Johanan L. Codona: 20051003:  Steward Observatory

theta = pi*(theta+eps);

I = real(2*(besselj(1,theta+eps)./(theta+eps)));
