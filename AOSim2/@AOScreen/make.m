function PS = make(PS,L0,fixLF)

% MAKE: Synthesize a new random phase screen.
%
% PS = make(PS,L0,fixLF)
%
% Renders a Kolmogorov AOPhaseScreen.
%
% INPUTS:
% PS: An AOPhaseScreen object.  You need to have set r0 (Fried's length)
%     before calling this routine.
% L0:    An Outer scale. (default is 30m)
% fixLF: If non-zero, patch up the low spatial-frequencies so that the
%        structure function behaves over a larger range.  (currently
%        doesn't do anything.)
%
% Written by: Johanan L. Codona, Steward Observatory CAAO
% June 25, 2002: First tests.
% July 12, 2002: Reworked for the AO classes.
% Feb. 25, 2003: Took care of high spatial frequency filtering.
% 20090413 JLCodona: AOSim2.

fprintf('Rendering %s %s.\n',class(PS),PS.name);

if(nargin<3)
	fixLF = false;
end

if(nargin<2)
	L0 = PS.L0;
end

K0 = 2*pi/L0;
K02 = K0^2;

% Synthetically construct the screen in K space.

PS = setDomain(PS,'k');

[kx,ky] = coords(PS);
[KX,KY] = meshgrid(kx,ky);

kmax = min(max(kx),max(ky));

% Start with complex gaussian random numbers.
zero(PS);
addNoise(PS,1.);

% Shape it so that it has the correct shape.
% PSD like k^5/3, sqrt(PSD) like k^5/6.

constant = (0.78)*2*pi; % Roger Angel's constant is in parenthesis.
rootarea = sqrt(prod(extent(PS)));

K2 = KX.^2 + KY.^2;
K = sqrt(K2);

clear KX KY

% NOTE: I have smoothed so that we should be able to shift by
% sub-pixels.  The filter limit is the max/2.5 to make it slightly
% better than Nyquist.
%   alpha = 4;
alpha = 11/3;
%   fudge = 1.3613;
fudge = 1.37;
% fudge = 1.32;
filt = ((K02+K2).^(-alpha/4)) .* smoothedge(kmax/2-K,kmax/50);
% filt = fudge * filt .* ((PS.r0.^(-5/6))) .* (rootarea * constant) ;
filt = filt .* ((PS.r0.^(-5/6))) .* (rootarea * constant) ;

filt(1,1) = 0.0;     % kill the DC term to make it zero-mean.

PS.grid_ = PS.grid_ .* filt;
tX(PS);

real(PS);

% estimate the current r0...
Nest = 5;
d1 = mean(mean(abs(PS.grid_ - circshift(PS.grid_,[0 Nest])).^2));
dx = PS.dx;
r00 = Nest*dx*(6.88/d1)^(3/5);

PS.grid_ = PS.grid_ * ((r00 / PS.r0)^(5/6) * PS.lambdaRef/2/pi);

PS.touched = false;

