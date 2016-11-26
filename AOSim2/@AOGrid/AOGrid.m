classdef AOGrid < matlab.mixin.Copyable  % formerly classdef AOGrid < handle
    % The main AOSim2 class.
    %
    % This is more important than you think...
    % "x" is dim 2
    % "y" is dim 1
    % lowercase is 1D
    % uppercase is 2D.
    % underscore means internal.
    %
    % The only exception to this is the pupil coordinates in AOAperture.
    %
    % Once you work with ij and xy long enough, you will want to scream.
    %
    % 20090407: JLCodona
    % 20090415 JLCodona.  Added fft method and new fftgrid_ usage model.
    % 20100514 JLCodona.  Added a coordinates caching scheme.
    % 20141120 JLCodona.  Made the objects deep-copyable.
    % Change comments are in the git log from now on.
    
    %% Properties
    properties(Constant=true, GetAccess='protected')
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
        FFTSize = [0 0]; % Recommended FFT Size. [0 0] is the AOGrid size.
        Offset = [0 0];  % This is a user Offset to move the object.
        % Offset differs from origin_ in that the latter is used to
        % optimize storage, while Offset is a user control.
        nanmap = 0;      % If when interpolating we get NaNs, this settable value is used to fill.
        verbosity = 0;  % set this to >0 for info and intermediate plots.
    
        interpolate_method = []; % quick or selected method.  Empty [] sets to qinterp2. Otherwise use the named method.
        seed = []; % Set this to empty for unseeded.  Set seed value for repeatable random numbers.
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
        
        % Coordinates cache
        % the actual cache...
        X_ = [];
        Y_ = [];
        % Cache update clues
        Xextremes = []; % This only gets set by COORDS.
        Yextremes = [];
        Nx_ = nan;
        Ny_ = nan;
        
    end
    
    %% Methods
    methods
        %% Constructors
        function obj = AOGrid(nxy)
            % function obj = AOGrid(nxy)
            % This is the basic AOGrid contructor.
            % The argument can be a pixel size (single number is square)
            % or it can be another AOGrid in which case it is used as a
            % template.  In that case no data is copied, just dimensions.
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
        
        function BLURB = describe(G)
            % BLURB = describe(G)
            % Returns a one-line description of the object.
            
            BLURB = sprintf('%s [%s: %dx%d]',G.name,class(G),G.size);
        end
        
        function G = touch(G)
            G.fftgrid_ = [];
        end
        
        function NXY = axisPixel(G,NXY)
            if(nargin>1)
                G.AXIS_PIXEL = NXY;
            end
            
            NXY = G.AXIS_PIXEL;
        end
        
        % NOTE: coords() is defined in an external file.
        function [X,Y] = COORDS(A,local)
            % [X,Y] = COORDS(A,local)
            % Returns the 2D coordinates of the grid A.
            if(nargin>1)
                [x,y] = coords(A,local);
            else
                [x,y] = coords(A);
            end
            
            % Okay.  Try to not have to compute this.
            
            % see if anything has changed since the caching...
            if(A.Nx_==length(x) && A.Ny_==length(y))
                
                %A.Xextremes==x([1 end])
                %A.Yextremes==y([1 end])
                
                % yeah, I know this is ugly.  What I really want is a call
                % like...
                % if(&&(Boolean_Array)) ...
                % If you are going to complain about elegance, provide a
                % better solution.
                if(prod(double([(A.Xextremes==x([1 end])) (A.Yextremes==y([1 end]))]))~=0)
                    % Looks like nothing has changed.  Use the cache!
                    X = A.X_;
                    Y = A.Y_;
                    %fprintf('<Using COORDS cached values.>\n');
                    return;
                else
                    % fprintf('<COORDS Cache failed values test. x(%g %g) y(%g %g)>\n',...
                    % A.Xextremes,x([1 end]),A.Yextremes,y([1 end]));
                end
            else
                % fprintf('<COORDS Cache failed length test. x(%g %g) y(%g %g)>\n',...
                % A.Nx_,length(x),A.Ny_,length(y));
            end
            % Try to speed this line up...
            %[X,Y] = meshgrid(x,y);
            % Note that coords returns row vectors
            
            %fprintf('DEBUG: Computing COORDS for %s <%s>\n',class(A),A.name);
            
            % Go ahead and compute it.
            X = ones(length(y),1)*x;
            Y = y'*ones(1,length(x));
            
            % save the cached coordinates.
            A.X_ = X;
            A.Y_ = Y;
            A.Nx_ = length(x);
            A.Ny_ = length(y);
            
            A.Xextremes = x([1 end]);
            A.Yextremes = y([1 end]);
        end
        
        function dom = domain(obj)
            dom = obj.domain_;
        end
        
        function c = axis(obj)
            c = obj.axis_;
        end
        
        function sz = size(obj,DIM)
            % SIZE = AOGRID.size([dim])
            % Return the size of the data grid in an AOGrid object.
            
            if(nargin<2)
                sz = size(obj.grid_);
            else
                sz = size(obj.grid_,DIM);
            end
        end
        
        function sz = resize(obj,varargin)
            % newsize = OBJ.resize(newsize);
            % Change the size of an AOGrid object's data grid.
            % (Assigning a new data grid does this automatically.)
            % Note that resize leaves the grid spacing unchanged.
            % Use OBJ.spacing(dx) to set the new spacing.
            
            switch length(varargin)
                case 0
                case 1
                    arg = varargin{1};
                    obj.axis_ = AOGrid.AXIS_FACE;
                    
                    if(isscalar(arg))
                        obj.grid_ = zeros([1 1]*arg);
                    else
                        obj.grid_ = zeros(arg(1:2));
                    end
                    
                    obj.AXIS_PIXEL = AOGrid.middlePixel(size(obj.grid_));
                    
                case 2
                    obj.axis_ = AOGrid.AXIS_FACE;
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
        
        function G = downsample(G,DEMAG)
            % AOGRID.downsample(DEMAG)
            % Decrease the sampling and size of an AOGrid.
            
            G.grid(downsampleCCD(G.grid,DEMAG,DEMAG)/DEMAG^2);
            G.spacing(G.spacing*DEMAG);
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
            % s = spacing(obj,varargin)
            
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
        
        %%
        function g = grid(obj,nugrid,mask) % TODO: make this fancier.
            %% function g = grid(obj,nugrid,[mask])
            % Returns the array inside of an AOGrid.
            % e.g.
            % grid = F.grid();
            %
            % You can set the grid to a new one by
            % F.grid(newarray);
            % The return value in the newgrid case is the object.
            %
            % Set the new value with a mask by
            % F.grid(newarray,which_pixels);
            %   >>> note that this does not coerce the shape.
            %
            % F.grid(newarray).show;
            % TODO: Fix the masked version of the call.
            
            if(nargin==1) % Just read out the value.
                g = obj.grid_;
            else % Set the value to the input...
                if(nargin==2) % the CLASSIC behavior...
                    nugrid = squeeze(nugrid);
                    nugrid = nugrid(:,:,1);
                    nugrid = squeeze(nugrid);
                    %if(size(obj)==size(nugrid))
                    if(prod(size(obj))==numel(nugrid))  %  Also allow vector assignments
                        obj.grid_(:) = double(nugrid(:)); % Does not force shape.
                        obj.touch;
                    else
                        obj.resize(size(nugrid));
                        obj.grid_ = double(nugrid);
                        obj.touch;
                        %error('different sized grid assignment not supported (yet)');
                    end
                else % This when a mask is specified...
                    obj.grid_(mask(:)) = nugrid(:);
                end
                g = obj; % Note that I return the object if setting the grid.
            end
        end
        
        function G = setMaskedGrid(G,NEWVALS,MASK)
            % G = setMaskedGrid(G,NEWVALS,MASK)
            % This should be part of how the G.grid(NEWVALS,MASK) function
            % works, but it is broken at the moment.
            
            G.grid_(MASK(:)) = NEWVALS(:);
            
        end
        
        %% 
        function G = justPhasor(G)
            % G = justPhase(G)
            % Replace the grid amplitudes with 1.  
            % Only phase remains.
            G.grid(exp(1i*angle(G.grid)));
        end
        
        %% This is not really useful anymore.  Center is the default alignment.
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
            % DEPRECATED.  This is part of AOSim(1).  It is essentially
            % meaningless now.
            % AOGrids are always in the x-domain and always "centered".
            % Centering means that the origin or reference point is in the
            % interior of the array, not at the (1,1) corner as is a fresh
            % FFT result.
            % I still keep this because it is in some of the AOSim(1) code.
            % *YOU* shouldn't need it.  
            
            % warning('AOGrid:DEPRECATED','Try to keep it in the x domain');
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
            % DEPRECATED.  This is part of AOSim(1).  It is essentially
            % meaningless now.
            % AOGrids are always in the x-domain and always "centered".
            % Centering means that the origin or reference point is in the
            % interior of the array, not at the (1,1) corner as is a fresh
            % FFT result.
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
        
        function g = checkFFTSize(g)
            % make sure FFTSize makes sense.
            
            if(isscalar(g.FFTSize)) % only one dim specified.
                g.FFTSize = g.FFTSize(1)*[1 1]; % Canonical 2D grid.
            end
            
            % Replace any 0 sizes with the actual grid size.
            SZ = g.size();
            DEFAULT = (g.FFTSize==0);
            g.FFTSize(DEFAULT) = SZ(DEFAULT);
            % there should no longer be any zeros in FFTSize.
        end
        
        function fgrid = fft(g,FFTSize)
            % Performs an FFT and returns the complex grid.
            %
            % This does not alter the grid_.
            % It caches the result in fftgrid_ or returns the cached result.
            %
            % Use kcoords and KCOORDS to get the spatial freq coordinates for fftgrid_.
            % Wavelength-dependent coords like theta are defined in AOField.
            %
            %NOTE BENE: The fftgrid_ is ALWAYS centered.
            %
            % If FFTSize is 0 then just use the native size of the AOGrid.
            %
            % Note: Use AOField.mkHalo to interpolate into the fftgrid_ array.
            % The coords there are unrelated to these, they are simply
            % angle (n.b. wavelength-dependent) for the returned image.
            % kcoords returns k-space coords for the extended (or
            % truncated) grid according to FFTSize.  Be careful.
            %
            % Conditions for (re)computing the FFT:
            % (a) if the fftgrid_ is empty.
            % (b) if fftgrid_ is not the same size as FFTSize.
            % (c) The AOGrid has been touched.
            %
            % otherwise just return the cached fftgrid_.
            
            %%
            if(nargin>1) % A new FFTSize was specified.
                g.FFTSize = FFTSize;
            end
            
            if(isscalar(g.FFTSize)) % only one dim specified.
                g.FFTSize = g.FFTSize(1)*[1 1]; % Canonical 2D grid.
            end
            
            g.checkFFTSize;
            
            if(sum(g.FFTSize~=size(g.fftgrid_))) % differs from the old size setting.  Clear cache.
                g.fftgrid_ = []; % Clear FFT cache
            end
            
            %% DO THE FFT
            if(isempty(g.fftgrid_)) % cache is cleared. RECOMPUTE FFT.
                %fprintf('DEBUG: Computing FFT and caching result.\n');
                
                g.FAXIS_PIXEL = AOGrid.middlePixel(g.FFTSize);
                
                % Each part does the fft and centers it.
                if(prod(double(g.FFTSize==g.size()))) % FFTSize IS the AOGrid size
                    g.fftgrid_ = circshift(...
                        fft2(circshift(g.grid_,1-g.AXIS_PIXEL)),...
                        g.FAXIS_PIXEL-1);
                else % FFTSize is NOT the AOGrid size
                    
                    if(sum(g.FFTSize>g.size())==2)
                        % the requested size is larger than the grid in both dims
                        g.fftgrid_ = padarray(g.grid(),g.FFTSize-g.size(),'post');
                        g.fftgrid_ = circshift(...
                            fft2(circshift(g.fftgrid_,1-g.AXIS_PIXEL)),...
                            g.FAXIS_PIXEL-1);
                    else %one or both dims are smaller than the grid
                        
                        % JLC 20101008: This code is probably not very efficient.
                        % I'm hoping it doesn't get used very often.
                        % Please feel free to contribute better code.

                        % TODO: (JLC 20150919) This is buggy. 
                        % It causes a x-y transpose in
                        % the fft.  I also wonder if there shouldn't be a
                        % warning if you set the FFT size to be smaller
                        % than the grid, since the only time I ever have
                        % FFTSize<size(g) is in error.  I usually use
                        % FFTSize to increase resolution when I start with
                        % a small grid, which is what I was intending.
                        % Using it to shrink the resolution is inneficient
                        % and probably not what the user meant to do.
                        fprintf('WARNING: the FFTSize is smaller than the grid.  Is this really what you intended to do?\n');
                        
                        DX = g.spacing;
                        
                        x1_ = 1:g.FFTSize(1);
                        x1_ = fftshift(x1_); % this takes care of special cases.
                        x1_ = (x1_-x1_(1))*DX(1);
                        
                        x2_ = 1:g.FFTSize(2);
                        x2_ = fftshift(x2_); % this takes care of special cases.
                        x2_ = (x2_-x2_(2))*DX(2);
                        
                        [X1_,X2_] = meshgrid(x1_,x2_);
                        
                        g.fftgrid_ = g.interpGrid(X2_,X1_); % backwards? forwards? BUG CHECK
                        g.fftgrid_(isnan(g.fftgrid_)) = 0; % inefficient zero-padding
                        
                        g.fftgrid_ = circshift(fft2(g.fftgrid_),g.FAXIS_PIXEL-1);
                    end 
                end                
            else
                fprintf('DEBUG: Returning cached FFT result.\n');
            end
            
            fgrid = g.fftgrid_;
        end
        
        function fgrid = ifft(g,FFTSize)
            % Performs an INVERSE FFT and returns the complex grid.
            %
            % (see AOGrid.fft() for more info.  This is basically the same
            % function.)
            %%
            if(nargin>1) % A new FFTSize was specified.
                g.FFTSize = FFTSize;
            end
            
            if(isscalar(g.FFTSize)) % only one dim specified.
                g.FFTSize = g.FFTSize(1)*[1 1]; % Canonical 2D grid.
            end
            
            g.checkFFTSize;
            
            if(sum(g.FFTSize~=size(g.fftgrid_))) % differs from the old size setting.  Clear cache.
                g.fftgrid_ = []; % Clear FFT cache
            end
            
            %% DO THE FFT
            if(isempty(g.fftgrid_)) % cache is cleared. RECOMPUTE FFT.
                %fprintf('DEBUG: Computing FFT and caching result.\n');
                
                g.FAXIS_PIXEL = AOGrid.middlePixel(g.FFTSize);
                
                % Each part does the fft and centers it.
                if(prod(double(g.FFTSize==g.size()))) % FFTSize IS the AOGrid size
                    g.fftgrid_ = circshift(...
                        ifft2(circshift(g.grid_,1-g.AXIS_PIXEL)),...
                        g.FAXIS_PIXEL-1);
                else % FFTSize is NOT the AOGrid size
                    
                    if(sum(g.FFTSize>g.size())==2)
                        % the requested size is larger than the grid in both dims
                        g.fftgrid_ = padarray(g.grid(),g.FFTSize-g.size(),'post');
                        g.fftgrid_ = circshift(...
                            ifft2(circshift(g.fftgrid_,1-g.AXIS_PIXEL)),...
                            g.FAXIS_PIXEL-1);
                    else %one or both dims are smaller than the grid
                        
                        % JLC 20101008: This code is probably not very efficient.
                        % I'm hoping it doesn't get used very often.
                        % Please feel free to contribute better code.
                        
                        DX = g.spacing;
                        
                        x1_ = 1:g.FFTSize(1);
                        x1_ = fftshift(x1_); % this takes care of special cases.
                        x1_ = (x1_-x1_(1))*DX(1);
                        
                        x2_ = 1:g.FFTSize(2);
                        x2_ = fftshift(x2_); % this takes care of special cases.
                        x2_ = (x2_-x2_(2))*DX(2);
                        
                        [X1_,X2_] = meshgrid(x1_,x2_);
                        
                        g.fftgrid_ = g.interpGrid(X2_,X1_); % backwards? forwards? BUG CHECK
                        g.fftgrid_(isnan(g.fftgrid_)) = 0; % inefficient zero-padding
                        
                        g.fftgrid_ = circshift(ifft2(g.fftgrid_),g.FAXIS_PIXEL-1);
                    end 
                end                
            else
                fprintf('DEBUG: Returning cached FFT result.\n');
            end
            
            fgrid = g.fftgrid_;
        end
        
        function DK = dk(g)
            % DK = dk(g):
            % This returns the fftgrid_ pixel sizes for the current FFTSize.
            % This is based on FFTSize, not the size of the currently allocated
            % AOGrid.fftgrid_.  Take note.
            
            g.checkFFTSize;
            
            DK = 2*pi./(g.FFTSize .* g.spacing_);
        end
        
        function [kx,ky] = kcoords(A)
            % This computes the coordinates of fftgrid_.
            % This is based on FFTSize, not the size of the currently allocated
            % AOGrid.fftgrid_.  Take note.
            % If the fftgrid_ hasn't been allocated, this routine will just
            % use FFTSize.  The pixel size is based on the FFTSize and the
            % pixel sizes in the space domain.  It's possible for the
            % fftgrid_ to not exist when this is computed, since it doesn't
            % get allocated until a request for an fft has been generated.
            % This should be okay, just make sure you call fft BEFORE you
            % expect there to be an answer in fftgrid_!

            A.checkFFTSize;
            SZ = A.FFTSize; 
            CEN = AOGrid.middlePixel(SZ);
            A.FAXIS_PIXEL = CEN;
            
            dk = A.dk;
            
            kx = ((1:SZ(1))-CEN(1))*dk(1);
            ky = ((1:SZ(2))-CEN(2))*dk(2);
        end
        
        function [KX,KY] = KCOORDS(A)
            [kx,ky] = kcoords(A);
            [KY,KX] = meshgrid(kx,ky);
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

        function A = rmMean(A)
            % A = A.rmMean;
            % Subtract the mean value from the grid.
            % Boring way: A.grid_ = A.grid_ - A.mean;
            % AOSim2 way: A - A.mean;
            A - A.mean;
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
        % dex = A.dex()
        % Returns the grid intensity values in log10 (decades).
            
            dex = log10(A.mag2);
        end
        
        function dex = ndex(A)
        % dex = A.ndex()
        % Returns the grid intensity values in log10 normalized to the max value.
            dex = A.mag2;
            dex = log10(dex/max(dex(:)));
        end
        
        function m = mean(A)
            % m = A.mean()
            % Returns the mean value of the grid.
            
            m = mean(A.grid_(:));
        end
        
        function v = var(A)
            % v = A.var()
            % Returns the variance of the grid.
            
            v = var(A.grid_(:));
        end
        
        function s = std(A)
            % sigma = A.std()
            % Returns the standard deviation of the grid.
            
            s = std(A.grid_(:));
        end
        
        function s = sigma(A)
            % sigma = A.sigma()
            % Returns the standard deviation of the grid.
            
            s = A.std;
        end
        
        function CEN = centroid(A)
        % CEN = A.centroid;
        % Compute the centroid of the AOGrid in physical units.
        
        [X,Y] = A.COORDS;
        M0 = mean(A.grid_(:));
        dG = A.grid - M0;
        M1x = mean(dG(:).*X(:))/M0;
        M1y = mean(dG(:).*Y(:))/M0;
        
        CEN = [M1x,M1y];
        end
        
        function s = sum(A)
        % s = sum(A)
        % Sum all of the pixels in the grid.
        
            s = sum(A.grid_(:));
        end
        
        function b = isX(G)
            % BOOLEAN = AOGRID.isX();
            % This method is left over from AOSim1 when grids could be in
            % configuration space or in Fourier space.  In AOSim2 all
            % AOGrids are assumed to be in configuration space.
            % This method is deprecated.
            b=(strcmp(G.domain_,AOGrid.DOMAIN_SPACE)==1);
        end
        
        function b = isK(G)
            % BOOLEAN = AOGRID.isK();
            % This method is left over from AOSim1 when grids could be in
            % configuration space (X) or in Fourier space (K).  
            % In AOSim2 all AOGrids are assumed to be in configuration space.
            % This method is deprecated.
            b=(strcmp(G.domain_,AOGrid.DOMAIN_FREQ)==1);
        end
        
        function b = isCentered(G)
            % BOOLEAN = AOGRID.isCentered();
            % This method is left over from AOSim1 when grids could be
            % centered on the middle of the array or on the [1,1] corner,
            % as in a Fourier transform result.
            % In AOSim2 all AOGrids are assumed to be centered.
            % This method is deprecated.
            b=(strcmp(G.axis_,AOGrid.AXIS_FACE)==1);
        end
        
        function yn = isCommensurate(a,b)
            % yn = isCommensurate(a,b)
            % Returns a boolean if two AOGrids have the same geometry.
            yn = false;
            
            if(strcmp(a.domain,b.domain)~=1)
                warning('AOGrid:DOMAIN','mismatched domains.');
                yn = false;
                return;
            end
            
            if(strcmp(a.axis,b.axis)~=1)
                a.center;
                b.center;
            end
            
            if( AOGrid.differ(size(a),size(b)) || ...
                    AOGrid.differ(spacing(a),spacing(b)) || ...
                    AOGrid.differ(origin(a),origin(b)) || ...
                    AOGrid.differ(a.Offset,b.Offset) )
                yn = false;
            else
                yn = true;
            end
        end
        
        %%
        function sgrid = subGrid(G,y,x)
            % function sgrid = subGrid(G,y,x)
            % Extract a smaller part of a grid.
            % x and y are 1-d lists of coordinates.
            % example>> SUBGRID = GRID.subGrid(10:20,30:40);
            %
            % Note that if x is omitted, x and y are assumed to be the same.
            
            if(nargin<3)
                sgrid = G.grid_(y,y);
            else
                sgrid = G.grid_(y,x);
            end
        end
        
        function G = setPixel(G,n1,n2,value)
        % G = setPixel(G,n1,n2,value)
            G.grid_(n1,n2) = value;
        end
        
        function G = bumpPixel(G,n1,n2,value)
            % G.bumpPixel(n1,n2,value)
            % Add to the (n1,n2) pixel value.
            
            G.grid_(n1,n2) = G.grid_(n1,n2) + value;
        end

        function G = setPixel1(G,n,value)
        % G = G.setPixel1(n,value)
        % Sets the nth pixel without coercing the shape.
            G.grid_(n) = value;
        end
        
        function G = bumpPixel1(G,n,value)
        % G = G.bumpPixel1(n,value)
        % Adds to the nth pixel without coercing the shape.
            G.grid_(n(:)) = G.grid_(n(:)) + value(:);
        end
        
        function G = stampCircle(G,CENTER,radius,value,INVERT)
            % G.stampCircle(CENTER,radius,[value=1],[INVERT=false])
            % Imprints a cirrcle on the grid.
            % CENTER is in real units (i.e. not pixels).
            % Note that the edge is smoothed by G.smooth.
           
            if(nargin<5)
                INVERT = false;
            end
            if(nargin<4)
                value = 1;
            end
            
            [X1,X2] = G.COORDS;
            R = sqrt((X1-CENTER(1)).^2+(X2-CENTER(2)).^2);
            if(INVERT)
                STAMP = value * smoothedge(R-radius,G.smooth);
            else
                STAMP = value * smoothedge(radius-R,G.smooth);
            end
            
            G.grid_ = G.grid_ .* STAMP;
        end
        
        function g = interpGrid(G,varargin)
            % g = interpGrid(G,varargin)
            % Returns a grid of G values interpolated to the spec'd coords.
            %
            % This command currently has two forms:
            %
            % g = G.interpGrid(XARRAY,YARRAY);
            % g = G.interpGrid(AOGRID);

            switch length(varargin)
                case 0
                    error('No coordinates to interpolate to?');
                case 1
                    arg1 = varargin{1};
                    if(isa(arg1,'AOGrid'))
                        [Xout,Yout] = COORDS(arg1);
                    else
                        error('I do not understand what you want to do with a single %s',class(arg1));
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
            
            if(isempty(G.interpolate_method))
                g = qinterp2(GX,GY,G.grid,Xout,Yout);
            else
                g = interp2(GX,GY,G.grid,Xout,Yout,G.interpolate_method);
            end
            % g(isnan(g)) = 0;
        end
        
        function G = plotC(G,gamma)
        % G = G.plotC(gamma)
        % Plot an AOGrid in complex.
            
            g = G.grid;
            if(isreal(g))
                % 				I = abs(g).^2;
                
                %minI = min(g(:));
                %maxI = max(g(:));
                
                [x,y] = coords(G);
                
                imagesc(x,y,g);
                axis square;
                axis xy;
                colorbar;
            else
                I = abs(g).^2;
                
                minI = min(I(:));
                maxI = max(I(:));
                
                if(maxI == 0)
                    maxI = 1;
                end
                
                if(minI/maxI > 0.99)
                    minI = 0;
                end
                
                
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
        
        function G = plotI(G,LIMS)
        % G = G.plotI([LIMS])
        % Plot an AOGrid intensity.
            
                [x,y] = G.coords;
                if(nargin<2)
                    imagesc(x,y,G.mag2);
                else
                    imagesc(x,y,G.mag2,LIMS);
                end

                axis square;
                axis xy;
                colorbar;
        end
        
        function G = plotDex(G,LIMS)
        % G = G.plotDex([LIMS])
        % Plot an AOGrid log10 intensity.
            
                [x,y] = G.coords;
                if(nargin<2)
                    imagesc(x,y,G.dex);
                else
                    imagesc(x,y,G.dex,LIMS);
                end

                axis square;
                axis xy;
                colorbar;
        end
        
        function delZ = del(G)
            % delZ = G.del()
            %
            % Compute the gradient of the grid.
            % Note that the first element of the answer is dZdx.
            % delZ = [dZ/dx,dZ/dy].
            
            dxy = G.spacing;
            [Zx,Zy] = gradient(G.grid_);
            delZ = [Zx/dxy(1),Zy/dxy(2)];
        end
        
        function RESULT = convolve(G,KERNEL)
            % RESULT = G.convolve(KERNEL);
        
            RESULT = conv2(G.grid_,KERNEL,'same');
        end
        
        function NUM = numel(g)
            % NUM = numel(g)
            % Returns the number of elements in the grid.
            
            NUM = numel(g.grid_); 
        end
        
        function G = add_rand(G,scale)
            % G = G.add_rand(scale)
            % Add uniform random real values to grid.
            
            if(nargin<2)
                scale = 1;
            end
            
            G + scale*(rand(G.size)-0.5);
        end
            
        function G = add_randn(G,scale)
            % G = G.add_rand(sigma)
            % Add Gaussian random real values to grid.
            
            if(nargin<2)
                scale = 1;
            end
            
            G + scale*randn(G.size);
        end
        
        function G = add_crandn(G,sigma)
            % G = G.add_rand(sigma)
            % Add Gaussian random complex values to grid.
            
            if(nargin<2)
                scale = 1;
            end
            
            G + (scale/sqrt(2))*(randn(G.size)+1i*randn(G.size));
        end
        
        %% Overloaded operators.
        function a = plus(a,b)
            % a = plus(a,b)
            % Add something to an AOGrid.
            
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
            % a = uminus(a)
            % Negate an AOGrid.
            
            a.grid_ = -a.grid_;
            touch(a);
        end
        
        function a = minus(a,b)
            % Subtract something from an AOGrid.
            
            if(~isa(a,'AOGrid'))
                error('operator call not in canonical form.');
            end
            if(isscalar(b))
                % fprintf('AOGrid scalar minus.\n');
                a.grid_ = a.grid_ - b;
                touch(a);
            elseif(isa(b,'AOGrid'))
                % fprintf('AOGrid-AOGrid.\n');
                if(isCommensurate(a,b))
                    a.grid_ = a.grid_ - b.grid;
                    touch(a);
                else
                    % fprintf('AOGrid-AOGrid: non-commensurate grids.\n');
                    [X,Y] = a.COORDS;
                    bg = b.interpGrid(X,Y);
                    bg(isnan(bg)) = 0;
                    a.grid_ = a.grid_ + bg;
                    touch(a);
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
                touch(a);
            elseif(isa(b,'AOGrid'))
                % fprintf('AOGrid .* AOGrid.\n');
                if(isCommensurate(a,b))
                    a.grid_ = a.grid_ .* b.grid;
                    touch(a);
                else
                    % fprintf('AOGrid .* AOGrid: non-commensurate grids.\n');
                    [X,Y] = a.COORDS;
                    bg = b.interpGrid(X,Y);
                    bg(isnan(bg)) = b.nanmap;
                    a.grid_ = a.grid_ .* bg;
                    touch(a);
                end
            end
        end
        
        function G = shiftPixels(G,pixels)
            % G = G.shiftPixels(PIXELS)
            % Perform a circular shift of the grid by the vector PIXELS.
            
            G.grid_ = circshift(G.grid_,pixels);
            G.touch;
        end
        
        function AOGRID = padBy(AOGRID,PADDING,PADVAL)
            % AOGRID = AOGRID.padBy(PADDING,[PADVAL])
            %
            % Pads the grid by PADDING pixels.
            
            if(numel(PADDING) < 2)
                PADDING(2) = PADDING(1);
            end
            
            if(nargin<3)
                PADVAL = 0;
            end
            
            AOGRID.grid(padarray(AOGRID.grid,PADDING,PADVAL,'both'));
            
        end
        
        function AOGRID = importFITS(AOGRID,FITSNAME,FRAME)
            % AOGRID = AOGRID.importFITS(FITSNAME,[FRAME_NUMBER])
            % Read in an AOGrid from a FITS file.
            % The grid will change size to accomodate the data.
            %
            % Note that I use MATLAB's FITS functions.
            
            if(nargin<3)
                FRAME = 1;
            end
            
            % I am going to keep it simple here by reading in a cube and
            % grabbing the part I need.  If this is too inefficient for
            % some reason, email me.
            
            CUBE = fitsread(FITSNAME);
            if(FRAME>size(CUBE,3))
                fprintf('Warning: The FITS file does not have that many image planes.\n');
                return;
            end
            
            AOGRID.grid(double(CUBE(:,:,FRAME)));
        end
        
        function G = multiplyRows(G,vector)
            % AOGrid G.multiplyRows(vector);
            % 
            % Multiply each row by the vector.
            % Only works if the vector is the same size as the grid row
            % size.
            
            vector = vector(:);
            LV = length(vector);
            SZ = G.size;

            if(LV ~= SZ(2))
                fprintf('AOGrid.multiplyRows Error: vector length (%d) and row size (%d) differ.\n',...
                    LV,SZ(2));
                return;
            else
                for n=1:SZ(1)
                    G.grid_(n,:) = G.grid_(n,:) .* vector';
                end
            end
        end
        
        function G = multiplyCols(G,vector)
            % AOGrid G.multiplyCols(vector);
            % 
            % Multiply each column by the vector.
            % Only works if the vector is the same size as the grid column size.
            
            vector = vector(:);
            LV = length(vector);
            SZ = G.size;

            if(LV ~= SZ(1))
                fprintf('AOGrid.multiplyCols Error: vector length (%d) and column size (%d) differ.\n',...
                    LV,SZ(1));
                return;
            else
                for n=1:SZ(2)
                    G.grid_(:,n) = G.grid_(:,n) .* vector;
                end
            end
        end
        
        %% Hilbert Space Operations 
        
        function G = rmModes(G,MODES)
            % G = rmModes(G,MODES)
            
            if(size(MODES,1) ~= G.numel)
                fprintf('AOGrid.rmModes: ERROR: MODES length (%d) is not the same as the grid (%d).\n',...
                    size(MODES,1),G.numel);
                return;
            end
            
            G.grid(G.grid_(:)-MODES*(MODES'*G.grid_(:)));
        end
        
        function G = shift0(G,SHIFT)
            % G = shift0(G,SHIFT)
            % Works like circshift on the grid but with zeros instead of
            % circular values.
            
            % Zero out the off-shift regions and then circshift.
            
            if(norm(SHIFT) == 0) % No shift required.
                return;
            end
            
            if(SHIFT(1)~=0)
                if(SHIFT(1)>0) % Lose the top
                    G.grid_(end-SHIFT(1)+1:end,:) = 0;
                else % Lose the bottom
                    G.grid_(1:1-SHIFT(1)-1,:) = 0;
                end
            end
            
            if(SHIFT(2)~=0)
                if(SHIFT(2)>0) % Lose the RHS
                    G.grid_(:,end-SHIFT(2)+1:end) = 0;
                else % Lose the LHS
                    G.grid_(:,1:1-SHIFT(2)-1) = 0;
                end
            end
            
            G.grid_ = circshift(G.grid_,SHIFT); 
        end
        
        function G = normalize(G)
            % G.normalize: Normalize the grid so its max value is 1+0i.
            
            [mx,loc] = max(abs(G.grid_(:))); % Note this returns the max complex value 
            if(mx~=0)
                G.grid_ = G.grid_/G.grid_(loc);
            end
        end
        
        function g_ = vec(G,MASK)
            % g_ = vec(G,[MASK])
            % G.vec returns the contents of the AOGrid grid as a column vector. 
            % G.vec(MASK) does the same but only selected values.
           
            if(nargin<2)
                g_ = G.grid_(:);
            else
                g_ = G.grid_(MASK(:));
            end
        end

        function N = norm(G)
            % N = norm(G)
            % Compute the L2 norm of the AOGrid data.
           
            N = norm(G.grid_);
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
            % STATIC: org = AOGrid.middlePixel(n)
            % Figure out which pixel would be the origin in an FFT of
            % length n.
            if(isa(n,'AOGrid'))
                n = n.size;
            end
            org = (n+2-mod(n,2))/2;
        end
        
        function G = multRows(G,vector)
            % G = multRows(G,vector);
            % 
            % Multiply each row by the vector.
            % Only works if the vector is the same size as the grid row
            % size.
            
            LV = numel(vector);
            SZ = size(G);

            if(LV ~= SZ(2))
                fprintf('AOGrid.multRows Error: vector length (%d) and row size (%d) differ.\n',...
                    LV,SZ(2));
                return;
            else
                parfor n=1:LV
                    G(:,n) = G(:,n) * vector(n); % Keep memory operations local.
                end
            end
        end

        function G = multCols(G,vector)
            % static G = AOGrid.multCols(G,vector);
            % 
            % Multiply each column by the vector.
            % Only works if the vector is the same size as the grid column size.
            
            vector = vector(:);
            LV = length(vector);
            SZ = size(G);

            if(LV ~= SZ(1))
                fprintf('AOGrid.multCols Error: vector length (%d) and column size (%d) differ.\n',...
                    LV,SZ(1));
                return;
            else
                parfor n=1:SZ(2)
                    G(:,n) = G(:,n) .* vector;
                end
            end
        end

        function M = RVmerge(M,V)
            % static MV = AOGrid.RVmerge(M,V)
            %
            % Merge the vector V into M in the combination M * (V .* X) = MV * X;
            
            M = bsxfun(@times,M,V.');
        end
        
        function M = LVmerge(M,V)
            % static VM = AOGrid.LVmerge(M,V)
            %
            % Merge the vector V into M in the combination V .* (M * X) = VM * X;
            
            M = bsxfun(@times,M,V);
        end
    end % static methods
end % of classdef

