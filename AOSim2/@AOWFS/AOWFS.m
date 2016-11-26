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
    properties(GetAccess='public',SetAccess='public')
        usePyr; % allows user to select pyramid (1) or SH (0)
    end
    properties(GetAccess = 'public', SetAccess = 'protected')
		ZernikeTip_ = {};        % Phase screen for x-tip Zernike mode (global tip)
		ZernikeTilt_ = {};       % Phase screen for y-tilt Zernike mode (global tilt)
	end
    
    methods
        % Constructor
        function WFS = AOWFS(APER,spacing,usePyramid,maskName,D)
            if nargin < 5, D = APER.estimateD(); end     % Should provide D for globalTipTilt method
            if nargin < 4  || isempty(maskName), maskName=0;end
            if nargin < 3  || isempty(usePyramid), usePyramid=0; end
            BBox = APER.BBox;
            SIDES = BBox(2,:)-BBox(1,:);
            Nacts = ceil(SIDES./spacing);
            WFS = WFS@AOGrid(Nacts);
            WFS = WFS@AODetector();
            WFS.usePyr = usePyramid;
            if maskName
                indpup = load(maskName);
                WFS.masked = wrapmask(indpup);
            else
                WFS.masked = false(size(WFS));
            end
            WFS.spacing([1 1].*spacing);
            %Isn't below unnecessary??  PMH commented out 091107
            %WFS.Offset=(mod(WFS.size,2)==0).*WFS.spacing/2; % even number subaps get bumped by half to center.
            WFS.masked = false(size(WFS));
            %WFS.masked = (APER.interpGrid(WFS)<0.75); % This might be better as a static parameter.
            WFS.centerOn(APER);
            WFS.Offset = WFS.Offset+(mod(WFS.size,2)==0).*WFS.spacing/2; % even number subaps get bumped by half to center.
            WFS.init(APER,0.2); % I just pick a masking threshold here.
            

			% Compute screens for projecting global aperture tip/tilt
			pupilMask = (APER.grid >= 0.5);
			pupilComplimentMask = (pupilMask == 0);
			% 1/sum(sum(pupilMask)) to normalize dot product by pupil area
			% second 2.0 since 1 unit of Noll tip/tilt mode is 2.0 units of OPD at r=1
			% 1/(D / 2.0) small angle approximation for slope angle in radians
			ttscale = 2.0 / ((D / 2.0) * sum(sum(pupilMask)));
			WFS.ZernikeTip_ = AOScreen(APER);                % need size and coordinates
			WFS.ZernikeTilt_ = AOScreen(APER);               % need size and coordinates
			WFS.ZernikeTip_.name = 'Tip';
			WFS.ZernikeTip_.zero;
			WFS.ZernikeTip_.addZernike(1, -1, ttscale, D);   % positive tip with increasing columns
			tempGrid = WFS.ZernikeTip_.grid;
			tempGrid(pupilComplimentMask) = 0;
			WFS.ZernikeTip_.grid(tempGrid);
			WFS.ZernikeTilt_.name = 'Tilt';
			WFS.ZernikeTilt_.zero;
			WFS.ZernikeTilt_.addZernike(1, 1, -ttscale, D);  % positive tilt with increasing rows
			tempGrid = WFS.ZernikeTilt_.grid;
			tempGrid(pupilComplimentMask) = 0;
			WFS.ZernikeTilt_.grid(tempGrid);
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
            if ~WFS.usePyr  %ie: if not using pyramid imported mask
                WFS.masked = false(WFS.size); 
                for nx=1:length(xw)
                xselect = abs(xf-xw(nx))<dxw;
                    for ny=1:length(yw)
                        yselect = abs(yf-yw(ny))<dyw;
                        f = field(yselect,xselect);
                        WFS.masked(ny,nx) = mean(abs(f(:)).^2)<thresh;
                    end
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
            if WFS.usePyr == 1
                WFS.sensePyramid(AOField(A).planewave * A);
            else
                WFS.sense(AOField(A).planewave * A); % lambda may be the default but it shouldn't matter.
            end
            
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
        
        function WFS = sensePyramid(WFS,F,useNoise)
            %Pyramid sensing code 
                       
            if(~isa(F,'AOField'))
                error('As a WFS, I am confused and frightened by %s inputs.',class(F));
            end
            
            if nargin < 3
                useNoise = 0;
            end
            %display(sprintf('sensePyramid useNoise is %1i\n',useNoise))
            
            field=F.grid;
            
            %Note: don't resample using fft2, because it will truncate
            
            %Saw aliasing and artifacts when transformed field to psf at
            %native resolution => pad array.  
            %Make sure the field object has an ODD number of elements,
            %otherwise fft2 will give an error.
                       
            %Pad field
            fieldPad=padarray(field,size(field)+1);
            
            if(mod(size(fieldPad,1), 2)==0 || mod(size(fieldPad,2),2)==0)
                display(' ')
                display('Field has an even number of elements.  FFT2 may produce an error.')
                display(' ')
                pause
            end
            
            
            psfAtTipPad = fftshift(fft2(fieldPad));
            
            %Match MMT orientation:
            %psfAtTipPad = rot90(psfAtTipPad,-2);
            %psfAtTipPad = flipud(psfAtTipPad);

            %check that padded psf looks OK
            %imagesc(log(abs(psfAtTipPad))); daspect([1 1 1]);
            %title('log psfAtTipPad')
            %pause
            
            
            Nelq = (size(psfAtTipPad)+1)/2;
            
            %Now, for the wobbling!
            ww =2;
            w(1,:) = [ ww  ww  -ww -ww];
            w(2,:) = [ ww -ww  -ww  ww];
    
            ee = size(psfAtTipPad);
    
        	V1 = 1 : (ww + 1);
            V2 = (ee(1) - ww) : ee(1) ;
            
            for i=1:length(w)
                
                clear psfAtTip
                psfAtTip = circshift( psfAtTipPad , w(:,i) );
                
                for j = 1:2
                    if w(j,i) < 0
                        d(j,:) = V2;
                    else
                        d(j,:) = V1;
                    end
                end
    
                psfAtTip(d(1,:),:)=0;
                psfAtTip(:,d(2,:))=0;
                
                %subplot(1,2,1);
                %imagesc(log(abs(psfAtTipPad))); daspect([1 1 1]); title('orig');
                %subplot(1,2,2);
                %imagesc(log(abs(psfAtTip))); daspect([1 1 1]); title('shifted');
                %pause
                
                %split PSF into quadrants
                %make sure each quadrant has an odd number of elements for FFT2?
                
                clear psfA psfB psfC psfD
                
                psfA = psfAtTip( 1:Nelq(1)  , 1:Nelq(2) );
                psfB = psfAtTip( 1:Nelq(1)  , Nelq(2):end );
                psfC = psfAtTip( Nelq(1):end , 1:Nelq(2) );
                psfD = psfAtTip( Nelq(1):end , Nelq(2):end );
                
                if(mod(size(psfA,1), 2)==0 || mod(size(psfA,2),2)==0)
                    display(' ')
                    display('Quadrant of field has an even number of elements.  FFT2 may produce an error.')
                    display(' ')
                    pause
                end

                %psfpyr = [psfA,psfB;psfC,psfD];
                %figure(2)
                %imagesc(log(abs(psfpyr))); daspect([1 1 1]);
                %drawnow;

                %transform back to pupil plane of WFS and downsample 
                %to correct number of subapps
                clear pupilA pupilB pupilC pupilD
                
                pupilAo = abs(ifft2(psfA)).^2 ; 
                pupilBo = abs(ifft2(psfB)).^2 ;
                pupilCo = abs(ifft2(psfC)).^2 ;
                pupilDo = abs(ifft2(psfD)).^2 ;

                %chop off extra pupil created during initial padding
                %chop off slightly less than 1/3 on each side b/c
                %mask has dark border - use scaling factor to adjust
                %sc = 31/29;  
                n= ceil(size(pupilAo)/3-.01*size(pupilAo,1));
                %n= floor(0.5 * ( size(pupilAo) - ceil( sc*size(pupilAo)/3 ) ));
                pupilA = pupilAo(n(1):end-n(1)+1 , n(2):end-n(2)+1 );
                pupilB = pupilBo(n(1):end-n(1)+1 , n(2):end-n(2)+1 );
                pupilC = pupilCo(n(1):end-n(1)+1 , n(2):end-n(2)+1 );
                pupilD = pupilDo(n(1):end-n(1)+1 , n(2):end-n(2)+1 );


                %clf; % Clear the Figure each time step
                %pupilpyr = [pupilA,pupilB;pupilC,pupilD];
                %figure(2)
                %imagesc(pupilpyr); daspect([1 1 1]);
                %drawnow;
                
                %interpolate to final WFS sampling
                xp=linspace(1,size(pupilA,2),length(WFS.masked));
                yp=linspace(1,size(pupilA,1),length(WFS.masked));
                %xp=linspace(1,size(pupilA,1),length(WFS.masked)-1);
                %yp=linspace(1,size(pupilA,2),length(WFS.masked)-1);
                [xxp,yyp]=meshgrid(xp,yp);
				%%%%%%%%%%%%%%%%%%%%%
				%FIX ME!!!
				%%%%%%%%%%%%%%%%%%%%
				
				%Don't use interp to downsample.  It will just find best fit point at the location you
				%request, not best fit value to area between points Need to upsample to
				%an integer ratio of the final size, then bin.  
				%Ex: Final Nsubap=25.  PupilA=91x91.  Upsample to 100x100 and do 4x4 binning (summing).  


                pupilAi(:,:,i) = interp2(pupilA, xxp,yyp);
                pupilBi(:,:,i) = interp2(pupilB, xxp,yyp);
                pupilCi(:,:,i) = interp2(pupilC, xxp,yyp);
                pupilDi(:,:,i) = interp2(pupilD, xxp,yyp);
                
                
            end
            
            pupilAint = sum(pupilAi,3);
            pupilBint = sum(pupilBi,3);
            pupilCint = sum(pupilCi,3);
            pupilDint = sum(pupilDi,3);
            
            %pupilint = [ pupilAint , pupilBint ; pupilCint, pupilDint ];
            %imagesc(pupilint); daspect([1 1 1]);
            %drawnow;
            
            
            %Then normalize the intensity in each pupil by the ratio of the 
            % (sum of all 4 pupils) / sum of original field intensity.
            %Do this so that photon noise is correct.

            fieldSum = sum(sum(abs(fieldPad)));
            pupilSumI = sum(sum(pupilAint+ pupilBint + pupilCint + pupilDint));
            pupilNorm = fieldSum / pupilSumI;

            pupilAint = pupilAint*pupilNorm;
            pupilBint = pupilBint*pupilNorm;
            pupilCint = pupilCint*pupilNorm;
            pupilDint = pupilDint*pupilNorm;
                
            if useNoise==1
                %Then worry about adding photon noise, read noise.
    			%TO DO:should maybe worry about making things integers...
                quantumeff = 0.9;
                rnoise = 3;

                pupilAint = pupilAint*1e-12*quantumeff;
                pupilBint = pupilBint*1e-12*quantumeff;
                pupilCint = pupilCint*1e-12*quantumeff;
                pupilDint = pupilDint*1e-12*quantumeff;
               % 
                pupilAint = 1e12*imnoise(pupilAint,'poisson');
                pupilBint = 1e12*imnoise(pupilBint,'poisson');
                pupilCint = 1e12*imnoise(pupilCint,'poisson');
                pupilDint = 1e12*imnoise(pupilDint,'poisson');

                % Read Noise
                pupilAint = pupilAint + randn(size(pupilAint))*rnoise ;
                pupilBint = pupilBint + randn(size(pupilBint))*rnoise ;
                pupilCint = pupilCint + randn(size(pupilCint))*rnoise ;
                pupilDint = pupilDint + randn(size(pupilDint))*rnoise ;
            end
            %calculate slopes for each subaperture
            %Account for pyramid LR&UD flip
            pupilSum = pupilAint+ pupilBint + pupilCint + pupilDint;
            xslopes = ((pupilAint + pupilCint)-(pupilBint + pupilDint)) ./ pupilSum;
            yslopes = ((pupilCint + pupilDint)-(pupilAint + pupilBint)) ./ pupilSum;
            
            
            WFS.grid_ = (xslopes+1i*yslopes).*(1-WFS.masked);  
            
            
        end
        
        function WFS = quiver(WFS,overplot)
            [XW,YW] = COORDS(WFS);
            CSLOPES = WFS.grid - WFS.bias;
            SLOPESy = real(CSLOPES); % Sigh.  Remember y is coord 1, 
            SLOPESx = imag(CSLOPES); % x is coord 2.
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

		function [tip, tilt] = globalTipTilt(WFS, ATMO)        % radians
			% Projection of ATMO onto Zernike modes gives tip/tilt in radians (use local copy)
			A = AOScreen(WFS.ZernikeTip_);    % need size and coordinates
 			A * ATMO;                         % dot product with Zernike mode (handles interpolation)
			tip = sum(sum(A.grid));
			% figure(6);
			% ATMO.show;
			% figure(4);
			% WFS.ZernikeTip_.show;
 			% figure(5);
 			% A.show;
			% drawnow;
			A = AOScreen(WFS.ZernikeTilt_);   % need size and coordinates
 			A * ATMO;                         % dot product with Zernike mode (handles interpolation)
			tilt = sum(sum(A.grid));
 			% input 'ATMO dot Tip: Press ENTER to Continue...'
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

