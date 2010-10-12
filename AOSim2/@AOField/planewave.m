function F = planewave(F,amplitude,THETA)

% PLANEWAVE: Set an AOField to a plane wave moving down the z axis.
%
% USAGE: F = planewave(F,[amplitude],[THETA])
%
% INPUTS:
% F: The AOField to work on.
% amplitude: the plane wave amplitude. (default=1)
% THETA: The vector angle of the source, in arcsecs.
%        (if this arg is missing, the planewave is on-axis).
%
% OUTPUTS:
% F: The modified AOField object.
%
% NOTE!!!!!  This leaves the field in the SPACE domain!
%
% SEE ALSO:
% AOGrid.
%
% Written by: Johanan L. Codona, Steward Observatory: CAAO
% July 12, 2002
% 20090413 JLCodona: AOSim2

% global env;

if(nargin < 2)
	amplitude = complex(1);
end

if(nargin < 3)
	ON_AXIS = true;
else
	ON_AXIS = false;
end

if(ON_AXIS)
	if(isX(F))
		F.constant(amplitude);
	else
		zero(F);
		F.setPixel(1,1,amplitude);
	end
	
else  % OFF AXIS
	
	F.tX;
	[X,Y] = COORDS(F);
	k = 2*pi/F.lambda;
	thx = THETA(1) / 206265;
	thy = THETA(2) / 206265;
	kappax = k*thx;
	kappay = k*thy;
	
	F.grid_ = amplitude .* exp(-1i.*(kappax.*X+kappay.*Y));
end

