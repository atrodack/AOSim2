function MTF = circularMTF(R,D)

% MTF = circularMTF(R,D)
% 
% A theoretical version of the circular aperture MTF.
% Useful for computing Strehl and making image deconvolution filters.
% 
% Note that the normalization is intended to give unity at R==0.
% 
% 20140409 JLCodona

rho = R/D;

% MTF = 0.5*(D)^2*real(acos(rho)-rho.*sqrt(1-rho.^2));
MTF = 2/pi*real(acos(rho)-rho.*sqrt(1-rho.^2));

