function PHOTONS = photonzCUBE(CUBE,N0,DATATYPE)
  
% photons = photonzCUBE(CUBE,N0,DATATYPE)
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
% June 7, 2013   Added CUBE processing.

PHOTONS = zeros(size(CUBE),DATATYPE);

for n=1:size(CUBE,3)
    PHOTONS(:,:,n) = photonz(CUBE(:,:,n),N0);
end




