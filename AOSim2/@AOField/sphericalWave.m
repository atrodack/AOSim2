function F = sphericalWave(F,amplitude,RANGE,CENTER,CUTOFF)

% AOField.sphericalWave(amplitude,RANGE,[CENTER],[CUTOFF])
%
% INPUTS:
% amplitude: the amplitude of the wave.
% RANGE: The distance to the center of the sphere in m.
% CENTER: The transverse coords of the wave center.
% Defaults to CENTER = [0 0].
% CUTOFF is the width of a Gaussian mask that limits the angular spectrum.
% It is measured in Fresnel scales and defaults to 25.
%
% 20150423 JLCodona


if(nargin < 4)
    CENTER = [0 0];
end

if(nargin < 5)
    CUTOFF = 25.;
end

[X,Y] = F.COORDS;
R2 = (X-CENTER(1)).^2 + (Y-CENTER(2)).^2;
FRESNEL = sqrt(F.lambda*RANGE);

% F.grid(amplitude .* exp((-1i.*F.k/2/RANGE)*R2));
% F.grid(amplitude .* exp((1i.*F.k/2/RANGE-1/(CUTOFF*FRESNEL)^2)*R2));
Rmax = FRESNEL^2/2/F.dx;
% Rmax = Rmax/3; % Be conservative.
F.grid(amplitude .* exp((1i.*F.k/2/RANGE)*R2).*smoothUP(1-R2/Rmax^2,1));

end
