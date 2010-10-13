function photons = photonz(IMG,N0)
  
% photons = photonz(IMG,N0);
%
% This function makes a photon noise version of a 
% field intensity.
%
% INPUTS: 
%   F: this
%   N0: How many photons are in the entire image.
% LOOKS: The number of dithered images used to build the image.
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


Sum0 = sum(IMG(:));

IMG = double(IMG)*(1e-12*N0/Sum0);

photons = 1e12*imnoise(IMG,'poisson');

return;



