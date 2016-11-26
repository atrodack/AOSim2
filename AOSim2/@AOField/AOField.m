classdef AOField < AOGrid
	% Segment = AOSegment();
	%
	% 20090412 JLCodona: AOSim2 version.
	
	%% Properties
	% Static Constants
	properties(Constant=true, SetAccess = 'private')
		VBAND = 0.5556e-6;
		RBAND = 0.7e-6;
		JBAND = 1.215e-6;
		HBAND = 1.654e-6;
		KBAND = 2.179e-6;
        LBAND = 3.8e-6;
		MBAND = 4.769e-6;
		NBAND = 10.472e-6;
        BrGamma = 2.165e-6; % Brackett Gamma
        HAlpha = 0.6563e-6; % Hydrogen Alpha 
        PaBeta =1.2818e-6;  % Paschen Beta
        OIII = 0.5007e-6;   % Oxygen III
		
		Rayleigh_LGS = 532e-9;
		Sodium_LGS = 589e-9;
		HeNe_Laser = 632.8e-9;
	end
	
	% Public properties
	properties(Access='public')
		lambda = AOField.HBAND;
		z = 0;
	end
	
	% Private
	properties(Access='private')
	end
	
	%% Methods
	methods
		function obj = AOField(ref)
			obj = obj@AOGrid(ref);
        end

        function [HALO,thx,thy] = mkHALO(F,FoV,dth) % angles in arcsecs.
            % Computes the FFT of part of the field and interpolates it to
            % the fftgrid_.  
            % There are coords for the interpolated grid and actual
            % coordinates for the FFT of the AOGrid.
            % KCOORDS returns the coordinates for 
            
            halo = F.fft; % do this first to initialize the buffers, etc.
            
            [kx,ky] = kcoords(F);
            
            % This is the size of the FFT pixels...
            dTH = F.dk/F.k*206265;
            thx_ = kx/F.k*206265;
            thy_ = ky/F.k*206265;
            
            if(nargin<3)
                dth = median(diff(thx_));
            end
            
            % this is get a slightly larger region to interpolate into.
            % No more extrapolation NaNs!
            SELx = abs(thx_)<=(FoV+2*dTH(1)); 
            SELy = abs(thy_)<=(FoV+2*dTH(2));
            thx_ = thx_(SELx);
            thy_ = thy_(SELy);
            [THX_,THY_] = meshgrid(thx_,thy_);
        
            Nth = ceil(FoV/dth);
            
            thx = (-Nth:Nth)*dth;
            thy = thx;
            [THX,THY] = meshgrid(thx,thy);
            
            halo = halo(SELx,SELy);
            
            HALO = qinterp2(THX_,THY_,halo,THX,THY);
            %HALO = interp2(THX_,THY_,halo,THX,THY,'cubic');
            
        end
        
        function [PSF,thx,thy] = mkPSF(F,FoV,dth,N0)
            % [PSF,thx,thy] = mkPSF(F,FoV,dth,N0)
            % N0 is optional.  If included then photonize the PSF.
            if(nargin<3)
                [PSF,thx,thy] = mkHALO(F,FoV);
            else
                [PSF,thx,thy] = mkHALO(F,FoV,dth);
            end
            PSF = abs(PSF).^2;
            PSF(isnan(PSF)) = 0;
            
            if(nargin>3)
            
                Sum0 = sum(PSF(:));
                PSF = double(PSF)*(1e-12*N0/Sum0);
                PSF = 1e12*imnoise(PSF,'poisson');
                
            end
        end
        
        function F = plotPSF(F,FoV,dexRange,dth)
        % F = plotPSF(F,FoV,dexRange,dth)
            if(nargin<3)
                dexRange = [-4 0];
            end

            if(nargin<4)
                [PSF,thx,thy] = F.mkPSF(FoV);
            else
                [PSF,thx,thy] = F.mkPSF(FoV,dth);
            end
            maxPSF = max(PSF(:));
            imagesc(thx,thy,log10(PSF/maxPSF),dexRange);
            axis square;
			axis xy;
        end
        
        function IMAGE = interferometer(F,ref)
            IMAGE = abs(F.grid+ref).^2;
        end

		% returns angles in radians.
		function [thx,thy] = thcoords(A)
			[kx,ky] = kcoords(A);
			k = 2*pi/A.lambda;
			thx = kx/k;
			thy = ky/k;
		end
		
		function [THX,THY] = THCOORDS(A)
			[KX,KY] = KCOORDS(A);
			k = 2*pi/A.lambda;
			THX = KX/k;
			THY = KY/k;
        end
		
		function a = mtimes(a,b)
		% function a = mtimes(a,b)
        % This is overloaded to compute smart "times" operations for AOFields.
        
			if(isa(b,'AOPhaseScreen'))
				if(isCommensurate(a,b))
					if(isPhase(b))
						% fprintf('DEBUG: AOField*AOPhaseScreen<phase>: commensurate grids.\n');
						a.grid_ = a.grid_ .* exp((1i*b.lambdaRef/a.lambda)*b.grid_);
					else
						% fprintf('DEBUG: AOField*AOPhaseScreen<phasor>: commensurate grids.\n');
						a.grid_ = a.grid_ .* b.grid_;
					end
				else
					[X,Y] = a.COORDS;
					
					if(isPhase(b))
						% fprintf('DEBUG: AOField*AOPhaseScreen<phase>: non-commensurate grids.\n');
						bg = exp((1i*b.lambdaRef/a.lambda)*interpGrid(b,X,Y));
					else
						% fprintf('DEBUG: AOField*AOPhaseScreen<phasor>: non-commensurate grids.\n');
						bg = interpGrid(b,X,Y);
					end
					
					bg(isnan(bg)) = 1;
					a.grid_ = a.grid_ .* bg;
				end
			end
			
			if(isa(b,'AOScreen'))
				if(b.mirror)
					MIRROR = 2;
				else
					MIRROR = 1;
                end
			
                if(b.conjugate)
                    FLIP = -1;
                else
                    FLIP = 1;
                end
                
				if(isCommensurate(a,b))
                    %fprintf('DEBUG: AOField*AOScreen: commensurate grids.\n');
					a.grid_ = a.grid_ .* exp((FLIP*MIRROR*2*pi*1i/a.lambda)*b.grid);
				else
					[X,Y] = a.COORDS;
 					%fprintf('DEBUG: AOField*AOScreen: non-commensurate grids.\n');
					bg = exp((FLIP*MIRROR*2*pi*1i/a.lambda) * b.interpGrid(X,Y));
					
					%bg(isnan(bg)) = 1;
					bg(isnan(bg)) = b.nanmap;
					a.grid_ = a.grid_ .* bg;
                end
                
			elseif(isa(b,'AOAtmo'))
% 				fprintf('DEBUG: AOField*AOAtmo at altitude %f.\n',a.z);
				% opl = b.OPL(a,a.z);
				a.grid_ = a.grid_ .* exp((2*pi*1i/a.lambda)*b.OPL(a,a.z));
			else
				a = mtimes@AOGrid(a,b);
            end
            
            a.touch;
        end
        
          function wavenumber = k(this)
            wavenumber = 2*pi/this.lambda;
        end    
        
        function F = setIntensity(F,star,band, bandWidth,integrationTime)
            if isa(star,'AOStar')
                if nargin < 5
                    integrationTime = 1;
                end
                bandIntFlux = bandWidth * star.GetNumPhotons(band,integrationTime);
            elseif isa(star,'double')
                bandIntFlux = star;
            else
                error('Confusing input to SetIntensity')
            end
            
            
            [xf,yf] = coords(F); % Fine grid of the field.

            fieldArea = abs(xf(end)-xf(1)) * abs(yf(end) - yf(1));
            
            numPixels = size(F.grid,1) * size(F.grid,2);
            
            pixelArea = fieldArea/numPixels;
            
            F.grid_ = F.grid_ * bandIntFlux * pixelArea;
        end
    end
end
