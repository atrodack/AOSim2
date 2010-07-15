function photons = photonize(F,N0,FOCALPLANE)
  
% photons = photonize(AOfield,N0,[FOCALPLANE]);
%
% This function makes a photon noise version of a field intensity.
%
% INPUTS: 
%   F: this
%   N0: How many photons are in the entire image.
%   FOCALPLANE: boolean for focal plane.
%
% OUTPUTS:
%   photons: The photon image.
%
% Written by: Johanan L. Codona, Steward Observatory: CAAO
% March 13, 2004
% March 19, 2004 Added flat field error (i.e. gain error).
% March 23, 2004 Reduced FF error by number of LOOKS.
% Dec 12, 2007   Rewrote this to work with a straight array and use 
%                MATLAB's imnoise function.
% 20090424 JLCodona.  Ported to AOSim2.

if(nargin>2 && FOCALPLANE)
    IMG = abs(F.fft).^2;
else
    IMG = F.mag2;
end




Sum0 = sum(IMG(:));
IMG = double(IMG)*(1e-12*N0/Sum0);
photons = 1e12*imnoise(IMG,'poisson');

return;

%% The following is the OLD CODE.
% 
% % photons = photonize(F,N0);
% %
% % This function makes a photon noise version of a field intensity.
% %
% % INPUTS: 
% %   F: this
% %   N0: How many photons are in the entire image.
% % LOOKS: The number of dithered images used to build the image.
% %
% % OUTPUTS:
% %   photons: The photon image.
% %
% % BUGS:
% %   I still don't have the statistics right for low-light pixels.
% %   I currently use a truncated version of the central limit
% %   theorem. If this is a problem, you should fix it!
% %
% % Written by: Johanan L. Codona, Steward Observatory: CAAO
% % March 13, 2004
% % March 19, 2004  Added flat field error (i.e. gain error).
% % March 23, 2004  Reduced FF error by number of LOOKS.
% 
% 
% if(nargin<3)
%   LOOKS = 1;
% end
% 
% FFERROR = 0.01/sqrt(LOOKS);
% 
% intensity = F.mag2;
% %max(max(intensity))
% %min(min(intensity))
% 
% intensity(isnan(intensity)) = 0;
% 
% Itotal = sum(intensity(:));
% intensity = intensity .* N0 ./ Itotal;
% 
% intensity = intensity .* (1+FFERROR*randn(size(intensity)));
% 
% Nmax = max(max(intensity));
% 
% photons = floor(max(0,intensity + sqrt(intensity) .* ...
% 		    randn(size(intensity))));
