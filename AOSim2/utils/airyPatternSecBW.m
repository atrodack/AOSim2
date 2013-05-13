function I = airyPatternSecBW(theta,obs,BW)

% function I = airyPatternSecBW(theta,obs,BW) : The Airy pattern intensity.
%
% INPUT:
%   theta: the set of angels to evaluate measured in lambda/D.
%
% OUTPUT:
%   I: The intensity at theta normalized to the peak.
%
% Johanan L. Codona: 20051003:  Steward Observatory
% 20080926: obs!  BW!!

thMax = max(theta(:));
stepMax = 0.1;
dScale = stepMax/thMax;

I = 0;
SCALE=demean(0:dScale:BW)+1;
for scale=SCALE
	I = I + airyPatternSecondary(scale*theta,obs);
end
I = I/length(SCALE);

