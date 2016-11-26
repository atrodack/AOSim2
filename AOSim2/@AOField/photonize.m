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
% 20090424 JLCodona.  Ported to AOSim2. Now an AOField method.

if(nargin>2 && FOCALPLANE)
    IMG = abs(F.fft).^2;
else
    IMG = F.mag2;
end

Sum0 = sum(IMG(:));
IMG = double(IMG)*(1e-12*N0/Sum0);
photons = 1e12*imnoise(IMG,'poisson');

return;
