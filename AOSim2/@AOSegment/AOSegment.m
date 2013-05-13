classdef AOSegment < AOGrid
    % Segment = AOSegment();
    %
    % 20090412 JLCodona: AOSim2 version.
    
    %% Properties
    % Public properties
    properties(Access='public')
        version = 1; % 1 is the classic pupil format.
        ignore = false; % don't include this segment in whatever.
        pupils; % for compatible use with version 1.
        smooth = 0.1;
        
        isMirror = true;  % false means piston has half the effect.
        piston = 0; % in m
        tiptilt = [0 0]; % in radians.
        lambdaRef = AOField.HBAND;
    end
    
    % Private
    properties(Access='private')
        pupilDef_ = {};
        touched = true;
    end
    
    %% Methods
    methods
        
        %% Constructors
        function obj = AOSegment(nxy)
            if(nargin==0)
                nxy = 1;
            end
            obj = obj@AOGrid(nxy);
			
			obj.nanmap = 0;  % make opaque beyond the pupil region.
			
        end
        
        %% Other methods
        function S = setPupilMode(S,md)
            fprintf('Changing pupilDef mode from %d to %d.\n',...
                S.version_,md);
            S.version_ = md;
        end
        
        function S = setPupil(S,pupil)
            fprintf('Setting pupil to new value.\n');
            S.pupilDef_ = {pupil};
        end
        
        function S = addPupil(S,pupil)
            fprintf('Appending another pupil defn.\n');
            S.pupilDef_ = {pupil};
		end
        
		function S = touch(S)
			touch@AOGrid(S);
			S.touched = true;
		end
		
        function g = grid(obj,nugrid) % TODO: make this fancier.
            if(obj.isMirror)
                OPL = 2*obj.piston;
            else
                OPL = obj.piston;
            end
            
            if(nargin==1)
                if(obj.tiptilt==0)
                    if(obj.piston==0)  % no TT,no piston
                        g = obj.grid_;
                    else  % no TT, just piston
                        g = exp(1i*2*pi*OPL/obj.lambdaRef) * obj.grid_;
                    end
                else  % piston and Tip tilt.
                    [X,Y] = COORDS(obj);
                    X0 = obj.Offset(1);
                    Y0 = obj.Offset(2);
                    g = exp((1i*2*pi/obj.lambdaRef)*...
                        (OPL+obj.tiptilt(1)*(X-X0) +...
                        obj.tiptilt(2)*(Y-Y0))) .* ...
                        obj.grid_;
                end
            else
                g = grid@AOGrid(obj,nugrid);
            end
        end
        
        function BB = BBox(S,local)
			% BBox is [x1min, x2min; x1max, x2max];
			% In other words...
			% [ ymin  xmin
			%   ymax  xmax ];
			% 
			% local is a flag that returns the BBox relative to the object
			% without the user Offset.  
            if(S.version == 1)
                if(isempty(S.pupils)) % for externally rendered segments.
                    [x,y] = S.coords;
                    filledy = (sum(S.grid,2) ~= 0);
                    filledx = (sum(S.grid,1) ~= 0);
                    
                    BB = [min(y(filledy)) min(x(filledx));
                        max(y(filledy)) max(x(filledx))];
                    return;
                end
                
                P = S.pupils;
				% Assume that the x y order in P is conventional (x,y).
				
                minX = min(P(:,1) - P(:,3)/2);
                maxX = max(P(:,1) + P(:,3)/2);
                
                maxY = max(P(:,2) + P(:,3)/2);
                minY = min(P(:,2) - P(:,3)/2);
                
                BB = [minY minX; maxY maxX];
                
                if(nargin>1 && local)
                    ORIGIN = [0 0];
                else
                    ORIGIN = S.Offset;
                end
                
                BB(:,1) = BB(:,1) + ORIGIN(1);
                BB(:,2) = BB(:,2) + ORIGIN(2);
            else
                error('I only understand AOSegment version 1 for now.');
            end
        end
        
        function S = plotBBox(S)
            BBOX = S.BBox;
            
            BB = [
                BBOX(1,2),BBOX(1,1);
                BBOX(1,2),BBOX(2,1);
                BBOX(2,2),BBOX(2,1);
                BBOX(2,2),BBOX(1,1);
                BBOX(1,2),BBOX(1,1);];
            
            hold on;
			AOSegment.plotBox(S.BBox);
            hold off;
        end
        
        function show(SEG) % overloading the AOGrid function.
            if(SEG.piston==0 && SEG.tiptilt(1)==0 && SEG.tiptilt(2)==0)
                show@AOGrid(SEG);
            else
                SEG.center;
                plotCAmpl(SEG.grid,1/2);
                title([class(SEG) ' ' SEG.name ': axis:' SEG.axis ' domain:' SEG.domain ],'FontSize',14);
                drawnow;
            end
        end
        
        function a = mtimes(a,b)
            a = mtimes@AOGrid(a,b);
        end
    end
    
    methods(Static=true)
        function plotBox(BBOX)
            BB = [
                BBOX(1,2),BBOX(1,1);
                BBOX(1,2),BBOX(2,1);
                BBOX(2,2),BBOX(2,1);
                BBOX(2,2),BBOX(1,1);
                BBOX(1,2),BBOX(1,1);];
            
            hold on;
            plot(BB(:,1),BB(:,2));
            hold off;
        end
    end
end
