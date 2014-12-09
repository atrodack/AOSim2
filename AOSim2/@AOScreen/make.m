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
% Merged rewrite into the git repository for regression testing.

fprintf('Rendering %s %s.\n',class(PS),PS.name);

% This enables a patch to correct for the structure function rolloff
% caused by using the intrinsically periodic FFT on a too-small grid.
if(nargin>2) 
	PS.fixLF = fixLF; % just setting this in the AOScreen should work too.
end

% This means a different outer scale was not defined.  Use the one in PS.
if(nargin<2)
	L0 = PS.L0;
end

% I like to do things in terms of spatial frequency.
% K0 is the spatial freq corresponding to the outer scale.
K0 = 2*pi/L0;
K02 = K0^2; % Just so I don't have to recompute it all the time.

% Synthetically construct the screen in K space.

% This uses old deprecated FFT features from the AOSim1 code. I probably 
% should change it.
PS = setDomain(PS,'k'); % sets to spatial freq domain without doing an FFT.
PS.axis_ = 'face';

% [kx,ky] = PS.kcoords;
% [KX,KY] = PS.KCOORDS;

SZ = PS.size;
CEN = AOGrid.middlePixel(SZ);

dk = PS.dk;

k1 = ((1:SZ(1))-CEN(1))*dk(1);
k2 = ((1:SZ(2))-CEN(2))*dk(2);

[KY,KX] = meshgrid(k2,k1);

kmax = min(k1(end),k2(end));

% Start with complex gaussian random numbers.
zero(PS);
addNoise(PS,1.);

% Shape it so that it has the correct shape.
% PSD like k.^(-5/3), sqrt(PSD) like k.^(-5/6).

% None of this matters since I am just going to normalize it at the end.

constant = (0.78)*2*pi; % Roger Angel's constant is in parenthesis.
rootarea = sqrt(prod(extent(PS)));

K2 = KX.^2 + KY.^2;
K = sqrt(K2);

clear KX KY

% NOTE: I have smoothed so that we should be able to shift by
% sub-pixels.  
% The filter limit is the max/2.5 to make it slightly better than Nyquist.

%   alpha = 4;
alpha = 11/3; % Kolmogorov
%   fudge = 1.3613;
% fudge = 1.37;
% fudge = 1.32;
filt = ((K02+K2).^(-alpha/4)); 
% filt = ((K02+K2).^(-alpha/4)) .* smoothedge(kmax/2.5-K,kmax/50);
% filt = ((K02+K2).^(-alpha/4)) .* smoothedge(kmax/2.5-K,kmax/50);
% filt = fudge * filt .* ((PS.r0.^(-5/6))) .* (rootarea * constant) ;
% filt = filt .* ((PS.r0.^(-5/6))) .* (rootarea * constant) ;

filt(1,1) = 0.0;     % kill the DC term to make it zero-mean.

PS.grid_ = PS.grid_ .* filt; % apply the filter.
tX(PS); % Transform back to the space domain.

% The screen is complex at this point. We just need the real part.

% real(PS);

SZ = PS.size;

if(PS.fixLF)
    Nlf = 8;
    %LFpart = interp2(imag(PS.grid_(1:round(SZ(1)/3),1:round(SZ(2)/3))),2,'cubic');
    %LFpart = interp2(imag(PS.grid_(1:round(SZ(1)/Nlf),1:round(SZ(2)/Nlf))),Nlf-1,'cubic');
    LFpart = imag(PS.grid_(1:round(SZ(1)/Nlf*1.1),1:round(SZ(2)/Nlf*1.1)));
    %LFpart = interp2(LFpart,Nlf-1);
    LFpart = interp2(LFpart,3);
%     whos LFpart
    LFpart = LFpart(1:SZ(1),1:SZ(2)); % trim to the size of the screen.
%     whos LFpart
    scaling = 0.8 * 8^(5/6);    
%     scaling = 5;    
    PS.grid_ = real(PS.grid_) + LFpart*scaling;
else
    PS.grid_ = real(PS.grid_);
end

% estimate the current r0...
Nest = round(PS.r0/PS.dx); % test near r0.
spacing = Nest*PS.dx;

% note that the screen is not necessarily periodic now.  
% Don't use circshift.

% d1 = mean(mean((PS.grid_ - circshift(PS.grid_,[0 Nest])).^2));
dZ = PS.grid_(:,1:end-Nest)-PS.grid_(:,1+Nest:end);
d1 = mean(mean(dZ.^2));
dx = PS.dx;
r00 = spacing*(6.88/d1)^(3/5);

PS.grid_ = PS.grid_ * ((r00 / PS.r0)^(5/6) * PS.lambdaRef/2/pi);

PS.touched = false;

