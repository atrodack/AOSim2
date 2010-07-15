classdef AOPhaseScreen < AOGrid
	% AOPHASESCREEN: The AOPhaseScreen class.
	%
	% USAGE:
	% PS = AOPhaseScreen(nx,ny,Cn2,thickness,altitude);
	%
	% INPUTS:
	% nx,ny:  Phase screen sizes. Empty inputs takes value from the
	% environment.
	%
	% Cn2:  The index structure constant for the layer. If empty or not
	%       specified, the value defaults to 3e-17 m^-2/3.
	%
	% thickness: thickness of the modeled turbulent layer (m) (Empty or not
	%       specified defaults to 100m).
	%
	% altitude: Height in meters of the layer.
	%
	% OUTPUTS:
	% PS: An AOPhaseScreen object.
	%
	% SEE ALSO:
	% AOGrid.
	%
	% Written by: Johanan L. Codona, Steward Observatory: CAAO
	% July 12, 2002
	
	
	% Static Constants
	properties(Constant=true, SetAccess = 'private')
		TYPE_PHASE = 'phase';
		TYPE_PHASOR = 'phasor';
	end
	
	% Public properties
	properties(Access='public', SetAccess='public')
		altitude = 0.;	% Default is on the ground.
		lambdaRef = AOField.VBAND;
	end
	
	% Private
	properties(SetAccess='private', GetAccess='public')
		touched = true;
		type = AOPhaseScreen.TYPE_PHASE;
		
		thickness;
		r0 = 0.15;
		Cn2;
	end
	
	%% Methods
	methods
		% Constructors
		function PS = AOPhaseScreen(nxy,r0)
			PS = PS@AOGrid(nxy);
			
			if(nargin>1)
				PS.r0 = r0;
			end
			
			PS.thickness = 100;
			PS.Cn2 = r0^(-5/3)/0.423/(2*pi/PS.lambdaRef)^2/PS.thickness;
		end
		
		% Utilities
		
		function b = isPhase(G)
            b=(strcmp(G.type,AOPhaseScreen.TYPE_PHASE)==1);
		end
		
		function a = expi(a)
			% EXPI: make into a phasor array via exp(i*grid).
			%
			% USAGE: PS = expi(PS);
			%  Note that this chages the type parameter from 'phase' to 'phasor'.  It
			%  is an error to try to expi a phasor field.
			%
			% Written by: Johanan L. Codona, Steward Observatory: CAAO
			% July 11, 2002
			
			if(~isPhase(a))
				return;
			end
			
			a.grid_ = exp(i*a.grid_);
			a.type = AOPhaseScreen.TYPE_PHASOR;
		end
		
		function PSI = phasor(PS)
			tX(PS);
			if(isPhase(PS))
				PSI = exp(1i*PS.grid_);
			else
				PSI = PS.grid_;
			end
		end
		
		function PHASE = phase(PS)
			tX(PS);
			if(isPhase(PS))
				PHASE = PS.grid_;
			else
				PHASE = angle(PS.grid_);
			end
		end
		
		function Rf = fresnel(a)
			Rf = sqrt(a.lambdaRef * a.altitude);
		end
		
% 		function a = mtimes(a,b)
% 			a = mtimes@AOGrid(a,b);
% 		end
		
		
	end % public methods.
end
