classdef AODM < AOScreen
	% AODM: Deformable Mirror class for the AOSim2 package.
	%
	% 20090419 JLCodona: First definitions.
	
	properties
		actuators=[]; % [x,y,z,segment_id,enabled]
		bconds = []; % boundary condition points.
	end
	
	methods
		% Constructors
		function DM = AODM(varargin)
			DM = DM@AOScreen(varargin);
			DM.mirror = false; %FIX ME!!!
			DM.defineBC(30,8);
		end
		
		% Operations
		function DM = setDM(DM,varargin)
			WF = varargin{1};
			
			switch class(WF)
				case 'AOAtmo'
					DM.grid_ = -WF.interpGrid(DM);
					
				case 'AOScreen'
					DM.grid_ = -WF.interpGrid(DM);
					% case 'AOField'
					% 		DM.grid_ = interpGrid(WF);
				otherwise
					error('AODM.match: Setting my surface using a %s frightens and confuses me.',class(WF));
			end
			DM.defineBC(30,8);
		end
		
		function DM = setActs(DM,varargin)
			WF = varargin{1};
			
			switch class(WF)
				case 'double'
					if(length(WF)~=1 && length(WF) ~= DM.nActs)
						error('actuator list length mismatch.');
					else
						DM.actuators(:,3) = WF;
					end
					% TODO: clip.
					touch(DM);
					
				case 'AOAtmo'
					DM.actuators(:,3) = -WF.interpGrid(DM.actuators(:,1),DM.actuators(:,2));
					% TODO: clip.
					touch(DM);
					
				case 'AOScreen'
					DM.actuators(:,3) = -WF.interpGrid(DM.actuators(:,1),DM.actuators(:,2));
					% TODO: clip.
					touch(DM);
					
					% case 'AOField'
					% ???? DM.grid_ = interpGrid(WF);
					% TODO: clip.
					% touch(DM);
					
				otherwise
					error('AODM.match: %ss confuse me.',class(WF));
			end
        end
        
        function DM = setSegmentActs(DM,varargin)
			WF = varargin{1};
            segNum=varargin{2};
			
			switch class(WF)
                %TODO: add double case
				%case 'double'
				%	if(length(WF)~=1 && length(WF) ~= DM.nActs)
				%		error('actuator list length mismatch.');
				%	else
				%		DM.actuators(:,3) = WF;
				%	end
				%	% TODO: clip.
				%	touch(DM);
					
                %TODO: add AOAtmo case
				%case 'AOAtmo'
				%	DM.actuators(:,3) = -WF.interpGrid(DM.actuators(:,1),DM.actuators(:,2));
				%	% TODO: clip.
				%	touch(DM);
					
				case 'AOScreen'
                    ACTS = getSegmentActs(DM,segNum);
					ACTS(:,3) = -WF.interpGrid(ACTS(:,1),ACTS(:,2));
                    DM.actuators(DM.actuators(:,4)==segNum,1:3)=ACTS;
					% TODO: clip.
					touch(DM);
					
					% case 'AOField'
					% ???? DM.grid_ = interpGrid(WF);
					% TODO: clip.
					% touch(DM);
					
				otherwise
					error('AODM.match: %ss confuse me.',class(WF));
			end
		end
		
		function ACTS = getSegmentActs(DM,segNum)
			ACTS = DM.actuators(DM.actuators(:,4)==segNum,1:3);
		end
		
		function TT = estimateSegmentTT(DM,SEGLIST)
			TT = zeros(length(SEGLIST),2);
			
			for segNum=SEGLIST
				ACTS = getSegmentActs(DM,segNum);
				X = ACTS(:,1); X = X - mean(X);
				Y = ACTS(:,2); Y = Y - mean(Y);
				Z = ACTS(:,3); Z = Z - mean(Z);
				clear ACTS;
				
				CCx = corrcoef(X,Z);
				CCx = CCx(1,2);
				CCy = corrcoef(Y,Z);
				CCy = CCy(1,2);
				
				COEFSx = polyfit(X,Z,1); SLOPEx = COEFSx(1);
				COEFSy = polyfit(Y,Z,1); SLOPEy = COEFSy(1);
				
				TT(segNum,:) = [SLOPEx*abs(CCx),SLOPEy*abs(CCy)];
			end
			
			% Take care of a special case...
			TT(isnan(TT)) = 0;
			
        end
		
        function DM = flatten(DM)
            DM.actuators(:,3) = 0;
            DM.touch;
        end

        function DM = poke(DM,n,val)
            DM.actuators(:,3) = 0;
            DM.actuators(n,3) = val;
            DM.touch;
        end

        function DM = addPoke(DM,n,val)
            DM.actuators(n,3) = DM.actuators(n,3) + val;
            DM.touch;
        end
        
		function DM = addFocus(DM,FL)
			[X,Y] = COORDS(DM);
			DM.grid_ = DM.grid + (X.^2 + Y.^2)/(2*FL);
        end
        
        function DM = addQuilt(DM,PtV,SPACING)
            %This simulates the quilting effect of discrete actuators
            %SPACING is the width of the Gaussian bump
            %PMH 090903
			[X,Y] = COORDS(DM);
            Xpos=DM.actuators(:,1);
            Ypos=DM.actuators(:,2);
            n=length(DM.actuators);
            for i=1:n
                R2=(X-Xpos(i)).^2+(Y-Ypos(i)).^2;
                DM.grid_ = DM.grid +PtV*exp(-(R2/(SPACING.^2)));
            end    
		end
		
		function DM = removeMean(DM)
			SELECT = DM.actuators(:,5)~=0;
			MEAN = mean(DM.actuators(SELECT,3));
			DM.actuators(SELECT,3) = DM.actuators(SELECT,3) - MEAN;
		end
		
		function DM = clip(DM,stroke)
            if(numel(stroke)==1)
                HIGH = stroke/2;
                LOW = -stroke/2;
            else
                HIGH = stroke(2);
                LOW =  stroke(1);
            end
            
			DM.actuators(:,3) = max(min(DM.actuators(:,3),HIGH),LOW);
        end
        
        function DM = addActs(DM,COORDS,segNum,CENTER)
			% This adds more actuators to the list.
			% NOTE: This is not the incremental update method.  See
			% bumpActs.
			% COORDS is a list of [x,y] actuator coords.
			% segNum is really just a label for plotting and sorting.  All
			% actuators are treated alike and together.
			
			if(nargin<3)
				segNum = nan;
			end
			
			if(nargin>3)
				if(isa(CENTER,'AOAperture'))
					segCenter = CENTER.segList{segNum}.Offset;
					
					COORDS(:,1) = COORDS(:,1) + segCenter(2);
					COORDS(:,2) = COORDS(:,2) + segCenter(1);
				else
					COORDS(:,1) = COORDS(:,1) + CENTER(2);
					COORDS(:,2) = COORDS(:,2) + CENTER(1);
				end
				
				% Always add in the overall offset of the DM...
				COORDS(:,1) = COORDS(:,1) + DM.Offset(2);
				COORDS(:,2) = COORDS(:,2) + DM.Offset(1);
			end
			
			for n=1:size(COORDS,1)
				nn = size(DM.actuators,1)+1;
				% (x,y,z,seg,enabled)
				DM.actuators(nn,1:2) = COORDS(n,:);
				DM.actuators(nn,3) = 0;
				DM.actuators(nn,4) = segNum;
				DM.actuators(nn,5) = true;
			end
			touch(DM);
		end
		
		function n = nActs(DM)
			n = size(DM.actuators,1);
		end
		
		function DM = defineBC(DM,radius,npoints)
			theta = (1:npoints)'/npoints*2*pi;
			DM.bconds = [cos(theta) sin(theta)] * radius;
			
			DM.bconds(:,1) = DM.bconds(:,1) + DM.Offset(2);
			DM.bconds(:,2) = DM.bconds(:,2) + DM.Offset(1);

			DM.bconds(:,3) = 0;
			touch(DM);
        end
		
        function DM = disableActuators(DM,THESE)
            DM.actuators(THESE,5) = 0;
        end

        function DM = enableActuators(DM,THESE)
            DM.actuators(THESE,5) = 1;
        end
        
		function DM = bumpActs(DM,VALUES)
			% This method adds the input VALUES to the current actuator
			% settings.
			if(length(VALUES) ~= DM.nActs)
				error('actuator list length mismatch.');
			end
			DM.actuators(:,3) = DM.actuators(:,3) + VALUES;
			
			% TODO: clip.
			touch(DM);
		end
		
		function DM = addRippleActs(DM,K,amp,phase)
			DM.actuators(:,3) = DM.actuators(:,3) + ...
				amp*cos(K(1)*DM.actuators(:,1) + K(2)*DM.actuators(:,2) + phase);
			% TODO: clip.
			touch(DM);
		end
		
		function LIST = listActuators(DM,seg)
			if(nargin<2)
				LIST = DM.actuators(:,1:3);
			else
				SELECT = DM.actuators(:,4)==seg && DM.actuators(:,5)~=0 ;
				LIST = DM.actuators(SELECT,1:3);
			end
		end
		
		function DM = plotActuators(DM,show_labels)
            if(nargin<2)
                show_labels = false;
            end
            hold on;            
			SELECT = DM.actuators(:,5)~=0;
			plot(DM.actuators(SELECT,1),DM.actuators(SELECT,2),'ro','MarkerSize',2);
			plot(DM.actuators(~SELECT,1),DM.actuators(~SELECT,2),'kx','MarkerSize',3);
			plot(DM.bconds(:,1),DM.bconds(:,2),'bs');
			daspect([1 1 1]);
            axis xy;
            if(show_labels)
                for n=1:DM.nActs
                    text(DM.actuators(n,1),DM.actuators(n,2),...
                        sprintf('%d',n),'FontSize',8);
                end
            end
            
            
            hold off;
		end
		
		function DM = plotRegions(DM)
			hold on;
			SELECT = DM.actuators(:,5)~=0;
			voronoi([DM.actuators(SELECT,1);DM.bconds(:,1)],...
				[DM.actuators(SELECT,2);DM.bconds(:,2)],'r');
			plot(DM.actuators(~SELECT,1),DM.actuators(~SELECT,2),'kx');
			plot(DM.bconds(:,1),DM.bconds(:,2),'bs');
			hold off;
			daspect([1 1 1]);
			axis xy;
		end
		
		function DM = render(DM)
			if(DM.touched)
				%fprintf('DEBUG: DM.rendering\n');
				[X,Y] = COORDS(DM);
				SEL = DM.actuators(:,5)~=0;
				DM.grid_ = ...
					griddata([DM.actuators(SEL,1);DM.bconds(:,1)],...
					[DM.actuators(SEL,2);DM.bconds(:,2)],...
					[DM.actuators(SEL,3);DM.bconds(:,3)],...
					X,Y,'cubic');
					% X,Y,'cubic');
				DM.grid_(isnan(DM.grid_)) = 0; % TODO: rethink extrapolation.
				DM.touched = false;
			else
				%fprintf('DEBUG: DM.render using CACHED\n');
			end
		end
		
		function a = uminus(a)
			if(a.nActs == 0)
				a.actuators(:,3) = -a.actuators(:,3);
			else
				uminus@AOScreen(a);
			end
			touch(a);
		end
		
		function g = grid(DM)
			if(DM.nActs == 0)
				g = DM.grid_;
				return
			end
			DM.render;
			g = DM.grid_;
		end
		
		function DM = show(DM,RANGE,MASK)
			% This method allows for a plot to be made with the common
			% requirement of an externally set range and an AOAperture
			% mask.
			
			G = real(DM.grid());
			[x,y] = coords(DM);
			
			if(nargin<2 || isempty(RANGE))
				Gmin = min(G(:))-eps;
                Gmax = max(G(:))+eps;
            else
                Gmin = RANGE(1)-eps;
                Gmax = RANGE(2)+eps;
            end
            
			if(nargin<3)
                MASK = 1;
            else
                MASK = real(MASK.grid());
            end

            
            
            imagesc(x,y,G.*MASK,[Gmin Gmax]);
            % plotCAmpl(AOG.grid_,1/4);
            axis square;
            axis xy;
            title([class(DM) ' ' DM.name ': axis:' DM.axis ' domain:' DM.domain ],...
                'FontSize',14);
            colorbar;
            drawnow;
            
            
			
		end
	end
end
