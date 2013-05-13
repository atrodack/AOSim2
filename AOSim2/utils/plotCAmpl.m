function plotCAmpl(CPLX,gamma,RANGE)

%function plotCAmpl(CPLX,gamma,[RANGE])
% note that the RANGE is pre-gamma.
%
% JLC!

CPLX = squeeze(CPLX);
% CPLX(isnan(CPLX)) = 0;

% I = CPLX.*conj(CPLX);
I = abs(CPLX).^2;
PHASE = angle(CPLX);

if(nargin<2)
    gamma = 1;
end

if(nargin==3)
    plotPhaseAmplitude(I.^gamma,PHASE,RANGE.^gamma);
else
    plotPhaseAmplitude(I.^gamma,PHASE);
end
