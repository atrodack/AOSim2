classdef AOAtmo < AOGrid
	% AOAtmo: AOAtmo class from AOSim2.
	%
	% This class holds AOScreens to define will hold the wavefront as meters of displacement, which makes
	% more sense.
	%
	% Written by: Johanan L. Codona, Steward Observatory: CAAO
	% Feb. 21, 2009
	% 20090418, JLCodona.  AOSim2.
	
	% properties
	properties(Access='public')
		layers = {};
		BEACON; % [x,y,z] a single point.
		time = 0;
		GEOMETRY = true;
		z = 0;
	end
	
	methods
		% Constructor
		function ATMO = AOAtmo(varargin)
			ATMO = ATMO@AOGrid(varargin);
		end
		
		% Operations
		function ATMO = addLayer(ATMO,screen,alt)
			n = length(ATMO.layers)+1;
			L = struct; % start building the layer.
			L.name = sprintf('Layer %d:%s',n,screen.name);
			if(nargin>2)
				screen.altitude = alt;
			end
			L.screen = screen;
			L.Wind = [0 0];
			
			L.ignore = false;
			
			ATMO.layers{end+1} = L;
			fprintf('AOAtmo now has %d layers.\n',length(ATMO.layers));
		end
		
		function ATMO = deleteLayer(ATMO,n)
			if(n>0 && n<=length(ATMO.layers))
				ATMO.layers(n) = [];
			else
				error('AOAtmo: cannot delete layer %d.',n);
			end
		end
		
		function ATMO = disp(ATMO)
			fprintf('%s %s: %d layers.\n',class(ATMO),...
				ATMO.name,ATMO.nLayers);
			
			for n=1:ATMO.nLayers
				disp(ATMO.layers{n});
			end
		end

		
		function n = nLayers(ATMO)
			n = length(ATMO.layers);
		end
		
		function ATMO = touch(ATMO)
			for n=1:ATMO.nLayers
				ATMO.layers{n}.screen.touch;
			end
        end
		
        function ATMO = make(ATMO)
        
            for n=1:ATMO.nLayers
                 ATMO.layers{n}.screen.make;
            end
        end
        
        
		function ATMO = setBeacon(ATMO,x,y,z)
			if(nargin==4)
				ATMO.BEACON = [ x y z ];
			else
				if(nargin==2)
					ATMO.BEACON = x;
				else
					error('requires separate x,y,z inputs or [x y z]');
				end
			end
		end
		
		function ATMO = useGeometry(ATMO,yn) % include geometry factors in results?
			ATMO.GEOMETRY = yn;
		end
		
		function ATMO = setObsTime(ATMO,t) % set observation time.
			ATMO.time = t;
		end
		
		function ATMO = setObsAltitude(ATMO,z) % set default observation altitude.
			ATMO.z = z;
		end
		
		function R = geomDistances(ATMO,X,Y,z) % Distances from the BEACON.
			if(nargin<4)
				z = 0;
			end
			% if(nargin<3)
			% error('not enough input coordinates.');
			% end
			
			if(isempty(ATMO.BEACON))
				error('The BEACON coordinates are undefined.');
			end
			
			if(isa(X,'AOGrid'))
				G = X;
				[X,Y] = COORDS(G);
				R = sqrt((X-ATMO.BEACON(1)).^2+(Y-ATMO.BEACON(2)).^2+(z-ATMO.BEACON(3)).^2);
			else
				R = sqrt((X-ATMO.BEACON(1)).^2+(Y-ATMO.BEACON(2)).^2+(z-ATMO.BEACON(3)).^2);
			end
		end
		
		function [X,Y] = scaleCone(ATMO,X,Y,z,znew,SOURCE) % Shrink to the SOURCE (defaults to BEACON).
			if(nargin<6)
				SOURCE = ATMO.BEACON;
			end
			if(nargin<5)
				error('not enough input coordinates.');
			end
			
			if(isempty(SOURCE))
				error('The BEACON coordinates are undefined.');
			end
			
			Lz = SOURCE(3) - z;
			Dz = znew - z;
			X = X + (SOURCE(1)-X)*(Dz/Lz);
			Y = Y + (SOURCE(2)-Y)*(Dz/Lz);
		end
		
		function g = grid(ATMO,nugrid)
			if(ATMO.nLayers == 0)
				g = grid@AOScreen(ATMO,nugrid);
				return;
			end
			
			[X,Y] = COORDS(ATMO);
			g = ATMO.OPL(X,Y,ATMO.z);
			
		end
		
		% ignore flag management
		function ATMO = ignoreAllLayers(ATMO)
			for n=1:ATMO.nLayers
				ATMO.layers{n}.ignore = true;
			end
		end
		
		function ATMO = includeAllLayers(ATMO)
			for n=1:ATMO.nLayers
				ATMO.layers{n}.ignore = false;
			end
		end
		
		% path integration functions
		function opl = OPL_(ATMO,X,Y,z,SOURCE)
			opl = zeros(size(X));
			
			% fprintf('DEBUG: OPL_ layers: ');
			for n=1:ATMO.nLayers
				% fprintf('%d ',n);
				if(~ATMO.layers{n}.ignore)
					zLayer = ATMO.layers{n}.screen.altitude;
					W = ATMO.layers{n}.Wind;
					t = ATMO.time;
					if(zLayer<=ATMO.BEACON(3))
						[XLayer,YLayer] = scaleCone(ATMO,X,Y,z,zLayer);
						if(ATMO.layers{n}.screen.touched)
							ATMO.layers{n}.screen.make;
						end
						opl_ = interpGrid(ATMO.layers{n}.screen,...
							XLayer-W(2)*t,YLayer-W(1)*t);
						opl_(isnan(opl_)) = 0;
						opl = opl + opl_;
					end
				end
			end
			%fprintf('\n');

			if(ATMO.GEOMETRY)
				opl = opl + geomDistances(ATMO,X,Y,z);
			end
		end
		
		function opl = OPL(ATMO,varargin)
			% OPL: valid arg lists... (Always uses ATMO BEACON.)
			%      Nothing. Uses internal phase screen.  Don't use this as
			%      it confuses things.
			%      some other AOGrid with coords.  Uses ATMO z.
			%      X,Y.  Uses ATMO z.
			%      X,Y,z.  Ignores ATMO z.
			if(nargin == 1)
				opl = ATMO.grid;
				return;
			end
			
			switch length(varargin)
				case 0
				case 1
					arg = varargin{1};
					
					if(isa(arg,'AOGrid'))
						[X,Y] = COORDS(arg);
						opl = ATMO.OPL_(X,Y,ATMO.z,ATMO.BEACON);
					else
						error('I dont know what to do with this argument yet.');
					end
					
				case 2
					arg1 = varargin{1};
					arg2 = varargin{2};
					
					if(isa(arg1,'AOGrid'))
						[X,Y] = COORDS(arg1);
						opl = ATMO.OPL_(X,Y,arg2,ATMO.BEACON);
					elseif(isfloat(arg1) && isfloat(arg2) )
						opl = ATMO.OPL_(arg1,arg2,ATMO.z,ATMO.BEACON);
					else
						
						error('I dont know what to do with this argument yet.');
					end
					
				case 3
					arg1 = varargin{1};
					arg2 = varargin{2};
					arg3 = varargin{3};
					
					if(isfloat(arg1) && isfloat(arg2) && isfloat(arg3) && isscalar(arg3))
						opl = ATMO.OPL_(arg1,arg2,arg3,ATMO.BEACON);
					else
						
						error('I dont know what to do with this argument yet.');
					end
					
				otherwise
					error('EWTF: too many arguments?  AOAtmo still BETA!');
			end
        end
        
        function r0 = totalFriedScale(ATMO,lambda)

            if(nargin<2)
                lambda = AOField.VBAND;
            end
            
            SLABS = 0;
            
            for n=1:length(ATMO.layers)
                SLABS = SLABS + ATMO.layers{n}.screen.thickness * ATMO.layers{n}.screen.Cn2;
            end
            
            r0 = (0.423*(2*pi/lambda)^2*SLABS)^(-3/5); % see e.g. Roddier or Fried or Tatarskii or ANYBODY!
        end
        
        
        
    end
end
