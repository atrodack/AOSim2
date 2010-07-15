classdef AOGrid < handle
    % The main AOSim2 class.
    %
	% This is more important than you think...
    % "x" is dim 2
    % "y" is dim 1
    % lowercase is 1D
    % uppercase is 2D.
    %
	% The only exception to this is the pupil coordinates in AOAperture.
	% 
    % 20090407: JLCodona
    % 20090415 JLCodona.  Added fft method and new fftgrid_ usage model.
    
    %% Properties
    properties(Constant=true, GetAccess='protected', SetAccess='private')
        SECS_PER_RADIAN = 206265;
        
        DOMAIN_SPACE = 'x';
        DOMAIN_FREQ = 'k';
        
        AXIS_CORNER = 'corner';  % FFT mode.
        AXIS_FACE = 'face'; % fftshift mode.
        AXIS_ARBITRARY = 'arb';  % use ORIGIN_ to find the center.
    end
    properties(GetAccess='public',SetAccess='public')
        name;
        defaultSize = [1 1]; % size of the default grid.
        FFTSize = 256*[1 1]; % recommended FFT size.
        Offset = [0 0];  % This is a user Offset to move the object.
		% Offset differs from origin_ in that the latter is used to
		% optimize storage, while Offset is a user control.
		nanmap = 0;      % If when interpolating we get NaNs, this settable value is used to fill.
        
        % Cache of coords to avoid meshgrid overload.
        x_ = [];
        y_ = [];
        X_ = [];
        Y_ = [];
    end
    
    %     properties(GetAccess = 'protected', SetAccess = 'protected')
    properties(GetAccess = 'public', SetAccess = 'protected')
        grid_;  % the actual array.
        fftgrid_ = []; % for cacheing the FFT.  Keep it in x! Doubles as ffttouch.
        FAXIS_PIXEL = [0 0]; % keep coordinated with fftgrid_.
        
        axis_ = AOGrid.AXIS_CORNER;  % what does (1,1) mean?
        domain_ = AOGrid.DOMAIN_SPACE;
        
        spacing_ = [0.04 0.04]; % default
        origin_ = [0 0]; % real world coords of the grid 'center'. This is private as opposed to Offset.
        AXIS_PIXEL = [1 1]; % keep coordinated with axis_.
    end
    
    %% Methods
    methods
        %% Constructors
        function obj = AOGrid(nxy)
            if(nargin==0)
                obj.grid_ = zeros(obj.defaultSize);
            else
                if(iscell(nxy))
                    nxy = nxy{1};
                end
                if(isa(nxy,'AOGrid'))
                    % obj = AOGrid.copyobj(nxy);
                    obj.name = ['copy of ' nxy.name];
                    obj.defaultSize = nxy.defaultSize;
                    obj.FFTSize = nxy.FFTSize;
                    obj.grid_ = nxy.grid;
                    obj.axis_ = nxy.axis_;
                    obj.domain_ = nxy.domain_;
                    obj.spacing_ = nxy.spacing_;
                    obj.origin_ = nxy.origin_;
                    obj.AXIS_PIXEL = nxy.axisPixel;
					obj.Offset = nxy.Offset;
                else
                    if(isscalar(nxy))
                        obj.grid_ = zeros(nxy,nxy);
                    else
                        obj.grid_ = zeros(nxy(1:2));
                    end
                end
            end
            
            obj.center; % default to center-faced.
        end
        
        %% Utilities
        
        function G = touch(G)
            G.fftgrid_ = [];
            x_ = [];
            y_ = [];
            X_ = [];
            Y_ = [];

        end
        
        function NXY = axisPixel(G,NXY)
            if(nargin>1)
                G.AXIS_PIXEL = NXY;
            end
            
            NXY = G.AXIS_PIXEL;
        end
        
        % NOTE: coords() is defined in an external file.
        function [X,Y] = COORDS(A,local)
            if(nargin>1)
                [x,y] = coords(A,local);
            else
                [x,y] = coords(A);
            end
            if(isempty(A.X_) || isempty(A.Y_))
                [A.X_,A.Y_] = meshgrid(x,y);
            end
            X = A.X_;
            Y = A.Y_;
        end
        
        function dom = domain(obj)
            dom = obj.domain_;
        end
        
        function c = axis(obj)
            c = obj.axis_;
        end
        
        function sz = size(obj)
            sz = size(obj.grid_);
		end
		
		function sz = resize(obj,varargin)
            switch length(varargin)
                case 0
                case 1
                    arg = varargin{1};
                    A.axis_ = AOGrid.AXIS_FACE;
                    
                    if(isscalar(arg))
                        obj.grid_ = zeros([1 1]*arg);
                    else
                        obj.grid_ = zeros(arg(1:2));
                    end
                    
                    obj.AXIS_PIXEL = AOGrid.middlePixel(size(obj.grid_));
                    
                case 2
                    A.axis_ = AOGrid.AXIS_FACE;
                    obj.grid_ = zeros([varargin{1} varargin{2}]);
                    obj.AXIS_PIXEL = AOGrid.middlePixel(size(obj.grid_));
                    
                otherwise
                    error('too many arguments');
            end
            
            sz = size(obj);
            obj.fftgrid_ = [];
		end
		
		function o = origin(obj,varargin)
            switch length(varargin)
                case 0
                case 1
                    arg = varargin{1};
                    if(isscalar(arg))
                        obj.origin_ = [1 1]*arg;
                    else
                        obj.origin_ = [arg(1) arg(2)];
                        obj.fftgrid_ = [];
                    end
                    
                case 2
                    obj.origin_ = [varargin{1} varargin{2}];
                    obj.fftgrid_ = [];
                    
                otherwise
                    error('too many arguments');
            end
            
            o = obj.origin_;
        end
        
        function G = setBBox(G,BBox,pad)
			% This makes a grid that is the size of the BBox and centered
			% on it.
			% BBox is [x1min, x2min; x1max, x2max];
			% In other words...
			% [ ymin  xmin
			%   ymax  xmax ];            
            midpoint = mean(BBox); 
            sides = BBox(2,:) - BBox(1,:);
            
            if(nargin>1)
                sides = sides + abs(pad);
            end
            
            pixels = ceil(sides ./ G.spacing);
            G.resize(pixels);
            G.origin(midpoint);
        end
        
        function n = nx(obj)
            n = size(obj.grid_,2);
		end
		
		function n = ny(obj)
            n = size(obj.grid_,1);
        end
        
        function val = dx(obj)
            val = obj.spacing_(2);
		end
		
		function val = dy(obj)
            val = obj.spacing_(1);
        end
        
        function s = spacing(obj,varargin)
            switch length(varargin)
                case 0
                case 1
                    arg = varargin{1};
                    if(isscalar(arg))
                        obj.spacing_ = [1 1]*arg;
                    else
                        obj.spacing_ = [arg(1) arg(2)];
                        obj.fftgrid_ = [];
                    end
                    
                case 2
                    obj.spacing_ = [varargin{1} varargin{2}];
                    obj.fftgrid_ = [];
                    
                otherwise
                    error('too many arguments');
            end
            
            s = obj.spacing_;
		end
		
		function D = extent(obj)
            D = size(obj.grid_) .* obj.spacing_;
		end
		
		function G = centerOn(G,C)
			if(~isa(C,'AOGrid'))
				error('I can only center an AOGrid object on another AOGrid object.');
			end
			
			G.Offset = C.Offset;
		end
		
		function DK = dk_(obj)
            DK = 2*pi./(size(obj.grid_) .* obj.spacing_);
        end
        
        function g = grid(obj,nugrid) % TODO: make this fancier.
            if(nargin==1)
                g = obj.grid_;
            else
                if(size(obj)==size(nugrid))
                    obj.grid_ = nugrid;
                    obj.fftgrid_ = [];
                else
                    error('different sized grid assignment not supported (yet)');
                end
                g = obj;
            end
        end
        
        function A = center(A)
            
            % CENTER: Move the grid to be face-centered.
            %
            % usage: center(AOGrid)
            %
            % For Fourier transform convenience, the array is usually centered
            % on the [1,1] pixel ('corner' centered.).  This is not handy for a number of
            % operations, including imaging.  Hence, center(AOGrid).
            %
            % SEE ALSO: corner, shift, shiftPixels.
            %
            % Written by: Johanan L. Codona, Steward Observatory CAAO
            % July 11, 2002
            % 20090407: JLCodona New version for AOSim2.
            
            if(strcmp(A.axis_,AOGrid.AXIS_FACE))
                return;
            end
            
            for dim=1:2
                A.grid_ = ifftshift(A.grid_,dim);
            end
            
            A.axis_ = AOGrid.AXIS_FACE;
            A.AXIS_PIXEL = AOGrid.middlePixel(size(A.grid_));
		end
		
		function A = corner(A)
            % CORNER: Center the grid on the 'corner' [1 1] pixel.
            %
            % usage: corner(AOGrid)
            %
            % For Fourier transform convenience, the array is usually centered
            % on the [1,1] pixel ('corner' centered.).  This is not handy for a number of
            % operations, including imaging.  Hence, center(AOGrid).
            %
            % SEE ALSO: center, shift, shiftPixels.
            %
            % Written by: Johanan L. Codona, Steward Observatory CAAO
            % July 11, 2002
            % 20090407: JLCodona New version for AOSim2.
            
            if(strcmp(A.axis_,AOGrid.AXIS_CORNER))
                return;
            end
            
            for dim=1:2
                A.grid_ = ifftshift(A.grid_,dim);
            end
            
            A.axis_ = AOGrid.AXIS_CORNER;
            A.AXIS_PIXEL = [1 1];
        end
        
        function g = tX(g)
            % TRANSFORMX: Fourier transform an AOGrid into x-space.
            if(strcmp(g.domain,AOGrid.DOMAIN_SPACE))
                return;
            end
            
            g = corner(g);
            g.grid_ = fft2(g.grid_)/prod(g.extent);
            g.domain_ = AOGrid.DOMAIN_SPACE;
            
            g = center(g);
            g.fftgrid_ = [];
		end
		
		function g = tK(g)
            warning('AOGrid:DEPRECATED','Try to keep it in the x domain');
            % TRANSFORMK: Fourier transform an AOGrid into k-space.
            if(strcmp(g.domain,AOGrid.DOMAIN_FREQ))
                return;
            end
            
            g = corner(g);
            g.grid_ = ifft2(g.grid_)*prod(g.extent);
            g.domain_ = AOGrid.DOMAIN_FREQ;
            
            g = center(g);
            g.fftgrid_ = [];
		end
        
		function g = checkFFTSize(g,FFTSize)
			% keepin' it square.  Complaints to JLC.
			FFTSize = max(FFTSize)*[1 1]; 
			
			Ngrid = max(size(g));
			Nfft = max(g.FFTSize);
			
			if(Ngrid > Nfft)
				FFTSize = 2^ceil(log2(Ngrid));
				g.FFTSize = FFTSize*[1 1];
			end
		end
		
        function fgrid = fft(g,FFTSize)
            % Performs FFT and returns complex grid.
            % This does not alter the grid_ if in x domain.
            % It caches the result or returns the cached result.
            % Use kcoords and KCOORDS to get the FFT dims.
            % Wavelength-dependent coords like theta are defined in
            % AOField.
            % NOTES: The fftgrid_ is ALWAYS centered.
            % 20090415 JLCodona.
            
			if(nargin>1)
				FFTSize = max(FFTSize)*[1 1]; % keepin' it square.  Complaints to JLC.

                if(g.FFTSize ~= FFTSize)
                    g.fftgrid_ = [];
				end
                
                if(length(FFTSize)>1)
                    g.FFTSize = FFTSize;
                else
                    g.FFTSize = FFTSize*[1 1];
				end
				
			g.checkFFTSize(FFTSize);
			
			else
				g.checkFFTSize(1);
			end

            if(isK(g))
                g.tX;
            end
            
			if(isempty(g.fftgrid_))
				%fprintf('DEBUG: Computing FFT and cacheing result.\n');
				
				g.center(); % prepare to pad.
				gpad = padarray(g.grid,g.FFTSize-g.size,'post');
				gpad = circshift(gpad,1-g.axisPixel);
				% Do our own padding to efficiently control phase tilt.
				g.fftgrid_ = ifft2(gpad)*(prod(g.spacing)*prod(g.FFTSize));
				% g.fftgrid_ = ifft2(g.grid(),g.FFTSize(1),g.FFTSize(2))...
				%    *(prod(g.spacing)*prod(g.FFTSize));
				
				for dim=1:2
					g.fftgrid_ = ifftshift(g.fftgrid_,dim);
				end
				g.FAXIS_PIXEL = AOGrid.middlePixel(size(g.fftgrid_));
				%             else
				%fprintf('DEBUG: Returning cached FFT result.\n');
			end
			
			fgrid = g.fftgrid_;
		end
		
		% This returns the value for the current FFTSize.
        function DK = dk(obj)
            DK = 2*pi./(obj.FFTSize .* obj.spacing_);
		end
		
		function [kx,ky] = kcoords(A)
            SZ = A.FFTSize;
            CEN = A.FAXIS_PIXEL;
            dk = A.dk;
            
            kx = ((1-CEN(1)):(SZ(1)-CEN(1)))*dk(1);
            ky = ((1-CEN(2)):(SZ(2)-CEN(2)))*dk(2);
		end
		
		function [KX,KY] = KCOORDS(A)
            [kx,ky] = kcoords(A);
            [KX,KY] = meshgrid(kx,ky);
		end
		
		function A = zero(A)
            % ZERO: Set an AOGrid to zero.
            A.grid_ = zeros(A.size);
            A.fftgrid_ = [];
		end
		
		function A = constant(A,val)
            % ZERO: Set an AOGrid to zero.
            A.grid_(:,:) = val;
            A.fftgrid_ = [];
		end
		
		function A = zeroNaNs(A)
            % ZERONANS: Set bad values to zero.
            A.grid_(isnan(A.grid_)) = 0;
            A.fftgrid_ = [];
        end
        
        function g = real(A)
            % if no output argument, modify the grid_.
            % Otherwise, don't change, just return the real part of the
            % grid_.
            if(nargout<1)
                A.grid_ = real(A.grid_);
                A.fftgrid_ = [];
            else
                g = real(A.grid_);
            end
		end
		
		function g = imag(A)
            % See real() method.
            if(nargout<1)
                A.grid_ = imag(A.grid_);
                A.fftgrid_ = [];
            else
                g = imag(A.grid_);
            end
		end
		
		function g = phase(A)
            if(nargout<1)
                A.grid_ = angle(A.grid_);
                A.fftgrid_ = [];
            else
                g = angle(A.grid_);
            end
		end
		
		function g = abs(A)
            if(nargout<1)
                A.grid_ = abs(A.grid_);
                A.fftgrid_ = [];
            else
                g = abs(A.grid_);
            end
		end
		
		function g = mag2(A)
            if(nargout<1)
                A.grid_ = (A.grid_).*conj(A.grid_);
                A.fftgrid_ = [];
            else
                % g = (A.grid_).*conj(A.grid_);
                g = abs(A.grid_).^2;
            end
        end
        
        function dex = dex(A)
            dex = log10(abs(A.grid_).^2);
		end
		
		function dex = ndex(A)
            dex = abs(A.grid_).^2;
            dex = log10(dex/max(dex(:)));
        end
        
        function m = mean(A)
            m = mean(A.grid_(:));
		end
		
		function v = var(A)
            v = var(A.grid_(:));
		end
		
		function s = std(A)
            s = std(A.grid_(:));
        end
        
        function b = isX(G)
            b=(strcmp(G.domain_,AOGrid.DOMAIN_SPACE)==1);
		end
		
		function b = isK(G)
            b=(strcmp(G.domain_,AOGrid.DOMAIN_FREQ)==1);
		end
		
		function b = isCentered(G)
            b=(strcmp(G.axis_,AOGrid.AXIS_FACE)==1);
		end
		
		function yn = isCommensurate(a,b)
            yn = false;
            
            if(strcmp(a.domain,b.domain)~=1)
                warning('AOGrid:DOMAIN','mismatched domains.');
				yn = false;
                return;
            end
            
            if(strcmp(a.axis,b.axis)~=1)
                %                 % TODO: Fix this yourself.
                %                 warning('AOGrid:ALIGN','mismatched axis alignment. try g.center.');
                %                 return;
                a.center;
                b.center;
            end
            
%             %if(size(a)~=size(b))
%             if(AOGrid.differ(size(a),size(b)))
%                 return;
%             end
%             
%             %if(spacing(a)~=spacing(b))
%             if(AOGrid.differ(spacing(a),spacing(b)))
%                 return;
%             end
%             
%             %if(origin(a)~=origin(b))
%             if(differ(origin(a),origin(b)))
%                 return;
%             end
%             
%             %if(a.Offset~=b.Offset)
%             if(differ(a.Offset,b.Offset))
%                 return;
%             end

                        %if(size(a)~=size(b))
            if( AOGrid.differ(size(a),size(b)) || ...
                AOGrid.differ(spacing(a),spacing(b)) || ...
                AOGrid.differ(origin(a),origin(b)) || ...
                AOGrid.differ(a.Offset,b.Offset) )
    			yn = false;
            else        
    			yn = true;
            end
        end
        
        function sgrid = subGrid(G,y,x)
            sgrid = G.grid_(y,x);
		end
		
		function G = setPixel(G,n1,n2,value)
            G.grid_(n1,n2) = value;
        end
        
        function g = interpGrid(G,varargin)
            switch length(varargin)
                case 0
                    error('No coordinates to interpolate to?');
                case 1
                    arg1 = varargin{1};
                    if(isa(arg1,'AOGrid'))
                        [Xout,Yout] = COORDS(arg1);
                    else
                        error('I do not understand what you want with a single %s',class(arg1));
                    end
                case 2
                    Xout = varargin{1};
                    Yout = varargin{2};
                otherwise
                    Xout = varargin{1};
                    Yout = varargin{2};
                    
            end
            
            G.center();
            [GX,GY] = COORDS(G);
            g = qinterp2(GX,GY,G.grid,Xout,Yout);
            %g = interp2(GX,GY,G.grid,Xout,Yout,'cubic');
            % g(isnan(g)) = 0;
        end
        
        function G = plotC(G,gamma)
            
            g = G.grid;
            if(isreal(g))
                % 				I = abs(g).^2;
                
                minI = min(g(:));
                maxI = max(g(:));
                
                [x,y] = coords(G);
                
                imagesc(x,y,g);
                axis square;
				axis xy;
                colorbar;
            else
                I = abs(g).^2;
                
                minI = min(I(:));
                maxI = max(I(:));
                
                I = (max(min(I,maxI),minI)-minI)/(maxI-minI);
                % I should now be 0<I<1
                
                if(nargin>1 && ~isempty(gamma))
                    I = I.^(1/gamma);
                end
                
                PHASE = angle(g) - 1.1;  % better alignment of colors.
                
                RGB(:,:,1) = I .* (1+cos(PHASE-pi/2))/2;
                RGB(:,:,2) = I .* (1+cos(2*PHASE-pi/3))/2;
                RGB(:,:,3) = I .* (1+cos(PHASE+pi/2))/2;
                
                [x,y] = coords(G);
                
                imagesc(x,y,RGB,[0 1]);
                axis square;
				axis xy;
            end
		end
		
		% Compute the gradient of the grid.
		% Note that the first element of the answer is dZdx.
		% delZ = [dZ/dx,dZ/dy].
		% Sorry for the confusion, but this is not simply compatable with
		% MATLAB. JLC. 20091005
		% TODO: Error handling!
		function delZ = del(G)
			dxy = G.spacing;
			%Z = G.grid;
			[Zx,Zy] = gradient(G.grid_);
			delZ = [Zx/dxy(1),Zy/dxy(2)];
		end
        
        %% Overloaded operators.
        function a = plus(a,b)
            if(~isa(a,'AOGrid'))
                error('operator call not in canonical form.');
            end
            if(isnumeric(b))
                %fprintf('AOGrid scalar sum.\n');
                a.grid_ = a.grid_ + b;
                a.fftgrid_ = [];
            elseif(isa(b,'AOGrid'))
                %fprintf('AOGrid+AOGrid.\n');
                if(isCommensurate(a,b))
                    a.grid_ = a.grid_ + b.grid;
                    a.fftgrid_ = [];
                else
                    %fprintf('AOGrid+AOGrid: non-commensurate grids.\n');
                    [X,Y] = a.COORDS;
                    bg = b.interpGrid(X,Y);
                    bg(isnan(bg)) = 0;
                    a.grid_ = a.grid_ + bg;
                    a.fftgrid_ = [];
                end
            end
        end
        
        function a = uminus(a)
            % fprintf('AOGrid unary minus.\n');
            a.grid_ = -a.grid_;
            touch(a);
        end
        
        function a = minus(a,b)
            if(~isa(a,'AOGrid'))
                error('operator call not in canonical form.');
            end
            if(isscalar(b))
                % fprintf('AOGrid scalar minus.\n');
                a.grid_ = a.grid_ - b;
                a.fftgrid_ = [];
            elseif(isa(b,'AOGrid'))
                % fprintf('AOGrid-AOGrid.\n');
                if(isCommensurate(a,b))
                    a.grid_ = a.grid_ - b.grid;
                    a.fftgrid_ = [];
                else
                    % fprintf('AOGrid-AOGrid: non-commensurate grids.\n');
                    [X,Y] = a.COORDS;
                    bg = b.interpGrid(X,Y);
                    bg(isnan(bg)) = 0;
                    a.grid_ = a.grid_ + bg;
                    a.fftgrid_ = [];
                end
            end
        end
        
        function a = mtimes(a,b)
            if(~isa(a,'AOGrid'))
                error('operator call is not in canonical form.');
			end
			
			if(isnumeric(b))
                % fprintf('AOGrid scalar product.\n');
                a.grid_ = a.grid_ * b;
                a.fftgrid_ = [];
            elseif(isa(b,'AOGrid'))
                % fprintf('AOGrid .* AOGrid.\n');
                if(isCommensurate(a,b))
                    a.grid_ = a.grid_ .* b.grid;
                    a.fftgrid_ = [];
                else
                    % fprintf('AOGrid .* AOGrid: non-commensurate grids.\n');
                    [X,Y] = a.COORDS;
                    bg = b.interpGrid(X,Y);
                    bg(isnan(bg)) = b.nanmap;
                    a.grid_ = a.grid_ .* bg;
                    a.fftgrid_ = [];
                end
            end
        end
        
        function G = shiftPixels(G,pixels)
            G.grid_ = circshift(G.grid_,pixels);
            G.fftgrid_ = [];
        end
    end % of methods
    
    %% static methods
    methods(Static=true)
        
        function yn = differ(mat1,mat2)
            % This is to work around an apparent bug in MATLAB.
            
            yn = prod(double(mat1(:)==mat2(:)))==0;
            
%             if(prod(double(mat1(:)==mat2(:))))
%                 yn = false;
%             else
%                 yn = true;
%             end
            
        end
        
        function copy = copyobj(obj)
            % Create a shallow copy of the calling object.
            copy = eval(class(obj));
            meta = eval(['?',class(obj)]);
            for p = 1: size(meta.Properties,1)
                pname = meta.Properties{p}.Name;
                try
                    eval(['copy.',pname,' = obj.',pname,';']);
                catch this
                    fprintf(['\nCould not copy ',pname,' ',this ,'.\n']);
                end
            end
        end
        
        function org = middlePixel(n)
            org = (n+2-mod(n,2))/2;
        end
        
        function array = normalize(array)
            % normalize: Normalize an array so its max value is unity.
            % USAGE: normed = AOGrid.normalize(array)
            %
            % Johanan L. Codona, Steward Observatory, CAAO
            % August 27, 2002 - Got tired of doing this manually!
            % 20090424 JLCodona: Added this as a static method.
            
            mx = max(array(:));
            if(mx~=0)
                array = array/mx;
            end
        end
        
    end % static methods
end % of classdef

