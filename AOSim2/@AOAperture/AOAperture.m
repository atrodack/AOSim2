classdef AOAperture < AOSegment
    % function AOA = AOAperture(nx,ny,pupils,smooth)
    % 20090412 JLCodona: AOSim2 version.
    
    % properties
    properties(Access='public')
        segList = {};
        combined; % AOSegment: doubles as AOGrid cache and touch.
    end
    
    methods
        function BB = BBox(S,n)
			% Return the BBox of segment n or the BBox hull if no n.
            BB = zeros(2);
            if(nargin<2)
                for n=1:length(S.segList)
                    BB_ = S.segList{n}.Segment.BBox();
                    OFFSET = S.segList{n}.Offset;
                    BB_(:,1) = BB_(:,1) + OFFSET(1);
                    BB_(:,2) = BB_(:,2) + OFFSET(2);
                    
                    BB(1,:) = min(BB(1,:),BB_(1,:));
                    BB(2,:) = max(BB(2,:),BB_(2,:));
                end
            else
                if(n<=length(S.segList))
                    BB = S.segList{n}.Segment.BBox();
                    OFFSET = S.segList{n}.Offset;
                    
                    BB(:,1) = BB(:,1) + OFFSET(1);
                    BB(:,2) = BB(:,2) + OFFSET(2);
                end
            end
        end
        
        function n = addSegment(A,Segment,Offset)
            n = length(A.segList)+1;
            s = struct;
            s.name = sprintf('Segment %d: %s',n,Segment.name);
            s.n = n;
            s.piston = 0;
            s.tiptilt = [0 0];
            s.Segment = Segment;
            if(nargin>2)
                s.Offset = Offset;
            else
                s.Offset = [0 0];
            end
            A.segList{end+1} = s;
            fprintf('Aperture now has %d segments.\n',length(A.segList));
            A.combined = [];
		end
		
		function A = removeSegment(A,n)
            A.segList(n) = [];
            fprintf('Aperture now has %d segments.\n',length(A.segList));
            A.combined = [];
		end
		
		function A = setTipTilts(A,TT) % TT is Nx2 list of radian tilts.
			if(size(TT,1) ~= length(A.segList))
				error('setTipTilts:Mismatched number of TipTilts and segments.');
			end
			if(size(TT,2) ~= 2)
                whos TT
				error('setTipTilts:TipTilts have to be a list of 2D vectors.');
			end
			
			for n=1:length(A.segList)
				A.segList{n}.tiptilt = TT(n,:);
			end
			
			touch(A);
		end
		
		function A = bumpTipTilts(A,TT) % TT is Nx2 list of radian tilts.
			if(size(TT,1) ~= length(A.segList))
				error('setTipTilts:Mismatched number of TipTilts and segments.');
			end
			
			for n=1:length(A.segList)
				A.segList{n}.tiptilt = A.segList{n}.tiptilt + TT(n,:);
			end
			
			touch(A);
		end
		
		function A = setPistons(A,PISTON) % PISTONS is a list of PISTON in m.
			if(size(PISTON,1) ~= length(A.segList))
				error('setPistons:Mismatched number of PISTONS (%d) and segments (%d).',...
                    size(PISTON,1),length(A.segList));
			end
			
			for n=1:length(A.segList)
				A.segList{n}.piston = PISTON(n);
			end
			
			touch(A);
		end
		
		function A = bumpPistons(A,PISTON) % PISTONS is a list of PISTON in m.
			if(size(PISTON,1) ~= length(A.segList))
				error('setPistons:Mismatched number of PISTONS (%d) and segments (%d).',...
                    size(PISTON,1),length(A.segList));
			end
			
			for n=1:length(A.segList)
				A.segList{n}.piston = A.segList{n}.piston + PISTON(n);
			end
			
			touch(A);
		end
		
		function A = setTipTilt(A,n,TT) % set ONE segment TipTilt.
			if(n > length(A.segList) || n<1)
				error('setTipTilt:segment number (%d) out of range.',n);
			end
			
			A.segList{n}.tiptilt = TT;
			touch(A);
		end
		
		function A = setPiston(A,n,PISTON) % Set ONE segment PISTON.
			if(n > length(A.segList))
				error('setPistons:segment number (%d) out of range.',n);
			end
			
			A.segList{n}.piston = PISTON;
			touch(A);
		end
		
		function S = plotBBoxes(S)
            hold on;
            for n=1:length(S.segList)
                S.plotBox(S.BBox(n));
            end
            hold off;
        end
        
        function [x,y] = segCoords(A,n)
            if(n<=length(A.segList))
                [x,y] = coords(A.segList{n}.Segment);
                OFFSET = A.segList{n}.Offset;
                x = x+OFFSET(1);
                y = y+OFFSET(2);
                % grid(A.segList{n}.Segment));
            else
                error('segCoords: There is no segment %d.',n);
            end
		end
		
		function g = segGrid(A,n)
            if(n<=length(A.segList))
                A.segList{n}.Segment.piston = A.segList{n}.piston;
                A.segList{n}.Segment.tiptilt = A.segList{n}.tiptilt;
                g = grid(A.segList{n}.Segment);
            else
                error('segGrid: There is no segment %d.',n);
            end
        end
        
        function A = make(A)
            A.render;
        end
        
        function BIGSEG = render(A)
            if(~isempty(A.combined))
                BIGSEG = A.combined;
				%fprintf('DEBUG: A.render: using cached value.\n');
                return;
			end
			
            %fprintf('DEBUG: A.render: Rebuilding complex aperture mask.\n');
            % Build the cached version for the current wavelength.
			
			A.cleanup();
            A.combined = AOSegment;
            A.combined.spacing(A.spacing);
            A.combined.setBBox(A.BBox,2*A.dx); % WARNING: padding is hardwired here.
            A.combined.zero;
            
			% Note: This has to be done every time since the Segments can be
			% shared.  Think of it as a special purpose shared calculator 
			% that needs to be tailored for each case and may have been
			% used by someone else since the last time.
            for n=1:length(A.segList)
                A.segList{n}.Segment.piston = A.segList{n}.piston;
                A.segList{n}.Segment.tiptilt = A.segList{n}.tiptilt;
                A.segList{n}.Segment.lambdaRef = A.lambdaRef;
				
				offset_save = A.segList{n}.Segment.Offset; % save in case important.
                A.segList{n}.Segment.Offset   = A.segList{n}.Offset;
                
				A.combined + A.segList{n}.Segment; % This actually does the rendering.
                % plotCAmpl(A.combined.grid_); sqar; % DEBUG!!!
				A.segList{n}.Segment.Offset = offset_save; % This just returns the zero from cleanup.
            end

            BIGSEG = A.combined;
			A.AXIS_PIXEL = A.combined.AXIS_PIXEL;
		end
        
        function g = grid(A,nugrid) % TODO: make this fancier.
            if(nargin>2)
                warning('AOAperture:SYNTAX','AOAperture.grid() I am ignoring the nugrid you passed.');
            end
            
            g = A.render.grid;
        end
        
        function A = touch(A)
            A.combined = [];
            A.fftgrid_ = [];
        end
        
        function sz = size(A)
            A.render;
            sz = A.combined.size;
		end

		function nx_ = nx(A)
			sz = size(A);
			nx_ = sz(2);
		end
		
		function ny_ = ny(A)
			sz = size(A);
			ny_ = sz(2);
		end
		
		function A = trueUp(A) % Align all segments.
            for n=1:length(A.segList)
                A.segList{n}.piston = 0;
                A.segList{n}.tiptilt = [0 0];
			end
			touch(A);
        end
        
        function NXY = axisPixel(G,NXY)
            if(isempty(G.combined))
               
                if(nargin>1)
                    G.combined.middlePixel(NXY);
                end
            end

            NXY = G.combined.axisPixel;
		end
		
		function [x,y] = coords(A,local)
            % Sorry, but this local business is confusing.  I don't think I
            % ever used it.  I should probably strip it out and deprecate
            % it. JLC 20090428.
            
            if(isempty(A.combined))
                A.make;
            end
            
            if(nargin>1)
                [x,y] = coords(A.combined,local);
            else
                [x,y] = coords(A.combined);
            end
			x = x + A.Offset(2);
			y = y + A.Offset(1);
		end
		
		function A = cleanup(A) 
			
			% Make sure all the individual defining segments are centered
			% relative to their own coordinate systems.  render() moves
			% them around.
			for n=1:length(A.segList)
				A.segList{n}.Segment.Offset = [0 0]; 
			end
			touch(A);
		end
		
		function show(A) % overloading the AOGrid function.
			A.center.plotC;
			% plotCAmpl(A.grid,1/2);
			title([class(A) ' ' A.name ': axis:' A.axis ' domain:' A.domain ],...
				'FontSize',14);
% 			drawnow;
		end
		
		function inside = isInside(A,POINTS,thresh)
			% This takes a list of [x1 x2] coordinates.
			if(nargin<3)
				thresh = 0.9;
			end
			
			vals = A.interpGrid(POINTS(:,1),POINTS(:,2));
			inside = vals>=thresh;
		end
	end
end
