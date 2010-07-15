classdef AOWFS < AOGrid & AODetector 
    % AOWFS Class
    % This is a simple non-physical WFS model.  Use it as a base
    % class for making better WFS models.
    %   20090421 JLCodona.  Part of AOSim2.
    
    properties
        masked = true;
        bias = 0; % n.b. I use complex numbers to represent the bias.
        qscale = 50; % for scaling the quiver plot arrows.
        NsubapFFT = nan;
    end
    
    methods
        % Constructor
        function WFS = AOWFS(APER,spacing)
            BBox = APER.BBox;
            SIDES = BBox(2,:)-BBox(1,:);
            Nacts = ceil(SIDES./spacing);
            WFS = WFS@AOGrid(Nacts);
            WFS = WFS@AODetector();
            WFS.spacing([1 1].*spacing);
            %Isn't below unnecessary??  PMH commented out 091107
            %WFS.Offset=(mod(WFS.size,2)==0).*WFS.spacing/2; % even number subaps get bumped by half to center.
            WFS.masked = false(size(WFS));
            %WFS.masked = (APER.interpGrid(WFS)<0.75); % This might be better as a static parameter.
            WFS.centerOn(APER);
            WFS.Offset = WFS.Offset+(mod(WFS.size,2)==0).*WFS.spacing/2; % even number subaps get bumped by half to center.
            WFS.init(APER,0.2); % I just pick a masking threshold here.
            

        end
        
        function N = nSubAps(WFS)
            N = sum(~WFS.masked(:));
		end
        
		function [ALL_SLOPES,MASK] = all_slopes(WFS)
            g = WFS.grid - WFS.bias;
            MASK = WFS.masked;
            ALL_SLOPES = [real(g(:));imag(g(:))];
		end
		
		function LIST = validSubapMap(WFS)
            VALID = ~WFS.masked;
			LIST = (1:length(VALID(:)))';
			LIST = LIST(VALID);
		end
		
        function WFS = init(WFS,A,thresh)
            % 20090820: Do a better job of setting masks and biases.
            
            [xw,yw] = coords(WFS);
            [xf,yf] = coords(A);
            
            field = A.grid;
            dxw = WFS.dx/2;
            dyw = WFS.dy/2;
            
            WFS.NsubapFFT = 2^ceil(log2(WFS.dx/A.dx));
            WFS.SetNumPixels(WFS.NsubapFFT);
            
            % Determine masking.
            WFS.masked = false(WFS.size);
            
            for nx=1:length(xw)
                xselect = abs(xf-xw(nx))<dxw;
                for ny=1:length(yw)
                    yselect = abs(yf-yw(ny))<dyw;
                    f = field(yselect,xselect);
                    WFS.masked(ny,nx) = mean(abs(f(:)).^2)<thresh;
                end
            end
            
            WFS = initBias(WFS,A);
        end
        
        %         function WFS = sense(WFS,F)
        %             % Compute the phase gradient over subaps.  I changed the
        %             % calculation to use medians instead of means due to phase
        %             % wrapping problems.  It is a lot better now. (JLC 20090422)
        %
        % 			% TODO: 20090819. Extend this to understand AOScreens. (MxH suggestion.)
        %
        % 			if(~isa(F,'AOField'))
        %                 error('As a WFS, I am confused and frightened by %s inputs.',class(F));
        % 			end
        %
        %             [xw,yw] = coords(WFS);
        %             [xf,yf] = coords(F);
        %             % [X,Y] = COORDS(WFS);
        %
        %             field = F.grid;
        % % 			kernel = chebwin(5);kernel=kernel*kernel';kernel=kernel/sum(kernel(:));
        % % 			field = conv2(field,kernel,'same');
        % % 			field(isnan(field)) = 0;
        %             dxw = WFS.dx/2;
        %             dyw = WFS.dy/2;
        %
        %             for nx=1:length(xw)
        %                 xselect = abs(xf-xw(nx))<=0.75*dxw;
        %                 for ny=1:length(yw)
        %                     if(WFS.masked(ny,nx))
        %                         WFS.grid_(ny,nx) = 0;
        %                     else
        %                         % The serious work happens here.
        %                         yselect = abs(yf-yw(ny))<=0.75*dyw;
        %                         f = field(yselect,xselect);
        %                         %f = angle(f/mean(f(:)));
        %                         f = angle(f);
        %                         dfdx = diff(f,1,2)/F.dx;
        %                         dfdy = diff(f,1,1)/F.dy;
        %                         %WFS.grid_(ny,nx) = mean(dfdx(:))+1i*mean(dfdy(:));
        %                         WFS.grid_(ny,nx) = median(dfdx(:))+1i*median(dfdy(:));
        % 						%imagesc(isnan(WFS.grid_));drawnow;
        %
        % 					end
        %                 end
        % 			end
        % 		end
        

		function WFS = initBias(WFS,A)
			% n.b. You may want to do A.trueUp before calling this if subaps can see across segments.
			WFS.sense(AOField(A).planewave * A); % lambda may be the default but it shouldn't matter.
            WFS.bias = WFS.grid;
		end

		
        function WFS = sense(WFS,F,useNoise)
            % Compute the subap centroid.
            % JLC 20090820: This is a more serious attempt at doing a Shack-Hartmann.
            % It still has a number of simplifications, but hopefully it
            % will have fewer runtime artifacts.
            % Remember to call init first.  The constructor should call it
            % automatically.
            
            % TODO: 20090819. Extend this to understand AOScreens. (MxH suggestion.)
            if(~isa(F,'AOField'))
                error('As a WFS, I am confused and frightened by %s inputs.',class(F));
            end
            
            if nargin < 3
                useNoise = false;
            end
            
            [xw,yw] = coords(WFS); % Coarse grid of the WFS subaps.
            [xf,yf] = coords(F); % Fine grid of the field.
            yf = yf';
            %dx = F.dx; % Size of the field pixels.  Assuming square. 
            
            field = F.grid;
            dxw = WFS.dx/2;
            dyw = WFS.dy/2;
            
            if ismethod(WFS,'SetNumPixels') && (~all(isnan(WFS.numPixels)) || numel(WFS.numPixels) ~= 1)
                N = WFS.numPixels(1);
            else
                N = WFS.NsubapFFT;
                WFS.SetNumPixels(N);
            end

            %N0 = N/2+1;
            %Was xx=-N/2+1:N/2;  This creates an offset -PMH
            xx = -N/2:N/2-1;
            
            for nx=1:length(xw)
                xselect = abs(xf-xw(nx))<dxw;
                for ny=1:length(yw)
                    if(WFS.masked(ny,nx))
                        WFS.grid_(ny,nx) = 0;
                    else
                        yselect = abs(yf-yw(ny))<dyw;
                        
                        f = field(yselect,xselect);
                        %psf = abs(fftshift2d(fft2(f,N,N))).^2 ;
                        psf = WFS.CreatePSF(f,useNoise);
                                                
                        norm = sum(psf(:));
                        C1 = mean(xx*psf)/norm;
                        C2 = mean(psf*xx')/norm;
                        WFS.grid_(ny,nx) = C1+1i*C2; % units are arbitrary.
                    end
                end
            end
        end
        
        function WFS = quiver(WFS,overplot)
            [XW,YW] = COORDS(WFS);
            CSLOPES = WFS.grid - WFS.bias;
            SLOPESx = real(CSLOPES);
            SLOPESy = imag(CSLOPES);
            if(nargin>1 && overplot~=0)
                hold on;
            end
            SLOPESx(isnan(SLOPESx)) = 0;
            quiver(XW,YW,10*SLOPESx/WFS.qscale,10*SLOPESy/WFS.qscale,1,'k');
            if(nargin>1 && overplot~=0)
                hold off;
            end
            daspect([1 1 1]);
        end
        
        function XY = subApCoords(WFS)
            [X,Y] = COORDS(WFS);
            X = X(~WFS.masked(:));
            Y = Y(~WFS.masked(:));
            XY = [X Y];
        end
        
        function SLOPES = slopes(WFS)
            g = WFS.grid - WFS.bias;
            valid = ~WFS.masked;
            valid = valid(:);
            SLOPES = [real(g(valid));imag(g(valid))];
        end
    end
    
    %% static methods
    methods(Static=true)
        
        function PISTONS = magicPistonSensor(F,A)
            N = length(A.segList);
            PISTONS = zeros(N,1);
            
            for n=1:N
                S = A.segList{n}.Segment;
                OFFSET_ = S.Offset; % keep.
                S.Offset = A.segList{n}.Offset;
                g = F.interpGrid(S);
                g(isnan(g)) = 0;
                g = g .* S.grid;
                S.Offset = OFFSET_; % restore.
                
                phasor_ = sum(g(:));
                phase = angle(phasor_);
                PISTONS(n) = phase/2/pi*F.lambda;
            end
            
        end
        
    end % static methods
    
end

