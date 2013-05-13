classdef AOReconstructor < handle
	%AOReconstructor class.
	%
	% 20090421: JLCodona.  UA.SO.CAAO.AOSim2
	
	properties
		A;
		DM;
		WFS;
		
		SLOPES;
		Scutoff = 10;
		
        % This is FYI only.  'fourier' or 'zernike' or 'diskharmonics' or 'kolmogorov'
		% ... or jlc_custom*, etc.
        TrainingMethod = 'fourier'; 
		% 	end
		%
		% 	properties(GetAccess = 'public', SetAccess = 'protected')
		RECONSTRUCTOR;
		
		ACTS;
		U,s,V; % SVD of SLOPES
		
		Nmodes = nan;
		
		D;
		OWD;
		step;
		lambda = AOField.RBAND;
		amplitude = AOField.RBAND/20;
        
        verbose = false; % print debugging info.
        
        d;       %Segment size for segment reconstructor
        CENTERS; %Segment centers for segment reconstructor
	end
	
	methods
		function RECON = AOReconstructor(A,DM,WFS)
			if(nargin<3)
				error('RECON = AOReconstructor(A,DM,WFS)');
			end
			
			RECON.A = A;
			RECON.DM = DM;
			RECON.WFS = WFS;
		end
		
		function RECON = program(RECON,D,OWD,step,lambda)
			% RECON = program(RECON,[D],[OWD],[step],[lambda])
			
			RECON.TrainingMethod = 'fourier';
			
			if(nargin<5)
				if(isempty(RECON.lambda))
					RECON.lambda = AOField.RBAND;
				end
			else
				RECON.lambda = lambda;
			end
			
			if(nargin<4)
				if(isempty(RECON.step))
					RECON.step = 0.5;
				end
			else
				RECON.step = step;
			end
			
			if(nargin<3)
				if(isempty(RECON.OWD))
					RECON.OWD = 6;
				end
			else
				RECON.OWD = OWD;
			end
			
			if(nargin<2)
				if(isempty(RECON.D))
					BBOX = RECON.A.BBox;
					RECON.D = mean(BBOX(2,:)-BBOX(1,:));
				end
			else
				RECON.D = D;
			end
			
			RECON.A.trueUp;
			RECON.WFS.initBias(RECON.A);
			
			dk = 2*pi/RECON.D * RECON.step;
			N = ceil(RECON.OWD/RECON.step);
			
            %Was -N:N, but this is redundant -PMH
			KSET = (0:N)*dk*RECON.step;
			KMAX = max(KSET);
			
			% We must do this for sin and cos.
			RECON.ACTS = zeros(RECON.DM.nActs,2*length(KSET)^2);
			RECON.SLOPES = zeros(2*RECON.WFS.nSubAps,2*length(KSET)^2);
			% 			disp(RECON);
			
			F = AOField(RECON.A);
			F.lambda = RECON.lambda;
			
			fprintf('AOReconstructor: Using ripples to feel out phase space...\n');
			
			n = 1;
			for kx=KSET
				for ky=KSET
					if(norm([kx ky])>KMAX)  % This forces the max spatial frequency probed to be the same in all directions.
						fprintf('x');
						continue;
					else
						fprintf('o');
					end
					RECON.DM.setActs(0).addRippleActs([kx ky],RECON.amplitude,0);
                    F.planewave * RECON.A * RECON.DM;
                    %Added RECON.A 20091107. Was line below previously
					%F.planewave * RECON.DM;
					acts = RECON.DM.actuators(:,3);
                    if RECON.WFS.usePyr == 1
                        RECON.WFS.sensePyramid(F);
                    else
                        RECON.WFS.sense(F);
                    end
					slopes=RECON.WFS.slopes;
					RECON.ACTS(:,n)   = acts;
					RECON.SLOPES(:,n) = slopes;
					n=n+1;
					
					RECON.DM.setActs(0).addRippleActs([kx ky],RECON.amplitude,pi/2);
					F.planewave * RECON.DM;
					acts = RECON.DM.actuators(:,3);
                    if RECON.WFS.usePyr == 1
                        RECON.WFS.sensePyramid(F);
                    else
                        RECON.WFS.sense(F);
                    end
					slopes=RECON.WFS.slopes;
					RECON.ACTS(:,n)   = acts;
					RECON.SLOPES(:,n) = slopes;
					n=n+1;
					
				end
				subplot(1,2,1);imagesc(RECON.ACTS);title('Actuators');
				subplot(1,2,2);imagesc(RECON.SLOPES);title('Slopes');
				drawnow;
				fprintf('\n');
			end
			fprintf('\n');
			
			% Trim off any unused space in the data arrays.
			RECON.ACTS(:,n:end)   = [];
			RECON.SLOPES(:,n:end) = [];
			
			fprintf('AOReconstructor: Using SVD to examine the slopes...\n');
			
			% NOTE and WARNING:
			% SVD pseudoinverses ONLY work from the left! 
			% Since we want to solve for a reconstructor on the left, we 
			% have to transpose first, do our thing, and transpose back.
			% Be advised. JLC
			
            RECON.SLOPES(isnan(RECON.SLOPES)) = 0;
            RECON.ACTS(isnan(RECON.ACTS)) = 0;
			
            [U,S,V] = svd(RECON.SLOPES','econ');
			RECON.U = U;
			RECON.s = diag(S);
			RECON.V = V;
			
			clear U S V
			
			RECON.rebuild;
		end
		
		function RECON = zprogram(RECON,D,Nmax,lambda)
			% RECON = zprogram(RECON,D,Nmax,lambda)
			
			if(nargin<5)
				if(isempty(RECON.lambda))
					RECON.lambda = AOField.RBAND;
				end
			else
				RECON.lambda = lambda;
			end
			
			RECON.TrainingMethod = 'zernike';
			RECON.OWD = Nmax;
			RECON.D = D;
			
			RECON.A.trueUp; 
			RECON.WFS.initBias(RECON.A);
			
			NZmodes = (Nmax+1)*(Nmax+2)/2;
			
			% We must do this for sin and cos.
			RECON.ACTS = zeros(RECON.DM.nActs,NZmodes);
			RECON.SLOPES = zeros(2*RECON.WFS.nSubAps,NZmodes);
			% 			disp(RECON);
			
			F = AOField(RECON.A);
			F.lambda = RECON.lambda;
			
			ABER = AOScreen(RECON.A);
			
			amp = RECON.amplitude;
			
			fprintf('AOReconstructor: Using ZERNIKES to feel out phase space to Order %d...\n',Nmax);
            fprintf('Number of modes being explored:%d\n',NZmodes);
			
			nmode = 1;
			
			for n=1:Nmax
                fprintf('\nOrder %d:',n);
                amp = 1*RECON.lambda/4/n;
				for m=-n:2:n
                    
                    weight = (n^2 + m^2).^(-5/6);
					%ABER.zero.addZernike(n,m,amp,RECON.D);
					ABER.zero.addZernike(n,m,weight*amp,RECON.D);
					RECON.DM.setActs(ABER);
					F.planewave * RECON.A * RECON.DM;
                    if RECON.WFS.usePyr == 1
                        RECON.WFS.sensePyramid(F);
                    else
                        RECON.WFS.sense(F);
                    end
					RECON.ACTS(:,nmode) = RECON.DM.actuators(:,3);
					RECON.SLOPES(:,nmode) = RECON.WFS.slopes;
					nmode = nmode + 1;
                    fprintf('%d ',m);
                    [x,y]=coords(F);
                    
                    if(RECON.verbose)
                        [x,y]=coords(F);
                        clf;
                        hold on;
                        subplot(1,2,1);
                        quiver(RECON.WFS,1);
                        %pause;
                        subplot(1,2,2);
                        imagesc(x,y,RECON.A.grid .* RECON.DM.grid); daspect([1 1 1]);
                        hold off;
                        drawnow;
                    end
				end
				%subplot(1,2,1);imagesc(RECON.ACTS);title('Actuators');
				%subplot(1,2,2);imagesc(RECON.SLOPES);title('Slopes');
                %drawnow;
% 				fprintf('\n');
			end
			fprintf('\n');
			
			fprintf('AOReconstructor: Using SVD to examine the slopes...\n');

            RECON.SLOPES(isnan(RECON.SLOPES)) = 0;
            RECON.ACTS(isnan(RECON.ACTS)) = 0;
            
			[UU,SS,VV] = svd(RECON.SLOPES','econ');
			RECON.U = UU;
			RECON.s = diag(SS);
			RECON.V = VV;
			
			clear UU SS VV
			
			RECON.rebuild;
		end
		
		function RECON = dhprogram(RECON,D,Nmax,verbose,lambda)
			% RECON = zprogram(RECON,D,Nmax,lambda)
			
			if(nargin<6)
				if(isempty(RECON.lambda))
					RECON.lambda = AOField.RBAND;
				end
			else
				RECON.lambda = lambda;
			end
			
			RECON.TrainingMethod = 'diskharmonics';
			RECON.OWD = Nmax;
			Mmax = round(pi*Nmax);
			RECON.D = D;
            RECON.verbose=verbose;
			
			NNmodes = Nmax*(2*Mmax+1);
			
			RECON.ACTS = zeros(RECON.DM.nActs,NNmodes);
			RECON.SLOPES = zeros(2*RECON.WFS.nSubAps,NNmodes);
			
			F = AOField(RECON.A);
			F.lambda = RECON.lambda;
			
			ABER = AOScreen(RECON.A);
			
			amp = RECON.amplitude; % scale these way down.
			
			if(amp>RECON.lambda/4/pi)
				fprintf('************************************************************************\n');
				fprintf('WARNING: The disk harmonic amplitudes use in training MAY be too large.\n');
				fprintf('************************************************************************\n');
			end
			
			fprintf('AOReconstructor: Using DISK HARMONICS to feel out phase space...\n');
			
			nmode = 1; % experiment counter.
			
			dh_init(); % Call Dr. Milton's code.
			
			RECON.WFS.qscale = 15;
			
			for n=1:Nmax % starting at 1  skips PISTON mode.
                fprintf('\nDisk Harmonic ORDER %d: ',n);
				for m=-Mmax:Mmax
					fprintf('%d ',m);
					
					ABER.zero.addDiskHarmonic(n,m,amp,RECON.D);
					RECON.DM.setActs(ABER);
					
					F.planewave * RECON.A * RECON.DM;
                    if RECON.WFS.usePyr == 1
                        RECON.WFS.sensePyramid(F);
                    else
                        RECON.WFS.sense(F);
                    end
                    
                    if(RECON.verbose)
                        [x,y]=coords(F);
                        clf;
                        hold on;
                        bigtitle(sprintf('DH mode(%d,%d)',n,m));
                        subplot(1,2,1);
                        quiver(RECON.WFS,1);
                        setFoV(D/2);
                        %pause;
                        subplot(1,2,2);
                        imagesc(x,y,RECON.A.grid .* RECON.DM.grid); 
                        setFoV(D/2);
                        daspect([1 1 1]);
                        hold off;
                        drawnow;
                    end
					
					RECON.ACTS(:,nmode) = RECON.DM.actuators(:,3);
					RECON.SLOPES(:,nmode) = RECON.WFS.slopes;
					nmode = nmode + 1;
                end
                
                if(RECON.verbose)
                    %subplot(1,2,1);imagesc(RECON.ACTS);bigtitle('Actuators');
                    %subplot(1,2,2);imagesc(RECON.SLOPES);bigtitle('Slopes');
                    %drawnow;
                    % fprintf('\n');
                end
                
			end
             fprintf('\n');
			
			fprintf('AOReconstructor: Using SVD to examine the slopes...\n');
			
			[U,S,V] = svd(RECON.SLOPES','econ');
			RECON.U = U;
			RECON.s = diag(S);
			RECON.V = V;
			
			clear U S V
			
			RECON.rebuild;
        end
                   
        function RECON = sprogram(RECON,D,Nmax,d,nmax,CENTERS,lambda)
            %This is a segment reconstructor.
            %It creates a few zernikes for the overall aperture
            %And then uses segment zernikes for the high order modes
			% RECON = sprogram(RECON,D,Nmax,lambda)
			%Segment-by-segment Interaction matrix
			if(nargin<8)
				if(isempty(RECON.lambda))
					RECON.lambda = AOField.RBAND;
				end
			else
				RECON.lambda = lambda;
			end
			
			RECON.TrainingMethod = 'segment-diskharmonic';
			RECON.OWD = Nmax;
			RECON.D = D;
            RECON.d= d;
            RECON.CENTERS=CENTERS;
			
			NZmodes = (Nmax+1)*(Nmax+2)/2;
			
			RECON.ACTS = zeros(RECON.DM.nActs,NZmodes);
			RECON.SLOPES = zeros(2*RECON.WFS.nSubAps,NZmodes);
			
			F = AOField(RECON.A);
			F.lambda = RECON.lambda;
			
			ABER = AOScreen(RECON.A);
			
			fprintf('AOReconstructor: Using ZERNIKES to feel out phase space to Order %2.0f...\n',Nmax);
            fprintf('Number of modes being explored:%3.0f\n',NZmodes);
			
			nmode = 1;
			
            %Full aperture zernikes
			for n=1:Nmax
                fprintf('\nOrder %d:',n);
                amp = 1*RECON.lambda/20;
				for m=-n:2:n
					ABER.zero.addZernike(n,m,amp,RECON.D);
					RECON.DM.setActs(ABER);
					F.planewave * RECON.A * RECON.DM;
                    if RECON.WFS.usePyr == 1
                        RECON.WFS.sensePyramid(F);
                    else
                        RECON.WFS.sense(F);
                    end
					RECON.ACTS(:,nmode) = RECON.DM.actuators(:,3);
					RECON.SLOPES(:,nmode) = RECON.WFS.slopes;
					nmode = nmode + 1;
                    fprintf('%d ',m);
                    [x,y]=coords(F);
                    clf;
                    hold on;
                    subplot(1,2,1); 
                    quiver(RECON.WFS,1);
                    %pause;
                    subplot(1,2,2); 
                    imagesc(x,y,RECON.A.grid .* RECON.DM.grid); daspect([1 1 1]);
                    hold off;
                    drawnow;
                    %pause;
                end
			end
			fprintf('\n');
            
            %Segment zernikes
            %Note the use of setsegmentActs to create segment aberrations
            segNum = length(RECON.A.segList);
            for a=1:segNum
               % for n=1:nmax
               %     ABER.zero;
               %     RECON.DM.setActs(ABER);
               %     fprintf('\nSegment %d, Order %d:',a,n);
               %     amp = 1*RECON.lambda/20;
               %     for m=-n:2:n
               %         ABER.zero.addZernike(n,m,amp,RECON.d,RECON.CENTERS(a,1),RECON.CENTERS(a,2));
               %         RECON.DM.setSegmentActs(ABER,a);
               %         F.planewave * RECON.A * RECON.DM;
               %         if RECON.WFS.usePyr == 1
               %            RECON.WFS.sensePyramid(F);
               %         else
               %            RECON.WFS.sense(F);
               %         end
               %         RECON.ACTS(:,nmode) = RECON.DM.actuators(:,3);
               %         RECON.SLOPES(:,nmode) = RECON.WFS.slopes;
               %         nmode = nmode + 1;
               %         fprintf('%d ',m);
               %         [x,y]=coords(F);
               %         clf;
               %         hold on;
               %         subplot(1,2,1); 
               %         quiver(RECON.WFS,1);
               %         %pause;
               %         subplot(1,2,2); 
               %         imagesc(x,y,RECON.A.grid .* RECON.DM.grid); daspect([1 1 1]);
               %         hold off;
               %         drawnow;
               %         %pause;
            %        end
               % end
    		%	fprintf('\n');
            
            Mmax = round(pi*nmax);
            amp = RECON.lambda/(8*pi); % really had to scale these down.
            
            %Segment disk harmonics
            %Note the use of setsegmentActs to create segment aberrations
			for n=1:nmax % starting at 1  skips PISTON mode.
                ABER.zero;
                RECON.DM.setActs(ABER);
                fprintf('\nTRAINING ON Disk Harmonic (%d)  Seg=%d.\n',n,a);
				for m=-Mmax:Mmax
					fprintf('%d ',m);
					RECON.DM.zero;
					ABER.zero.addDiskHarmonic(n,m,amp,RECON.d,RECON.CENTERS(a,2),RECON.CENTERS(a,1));
                    RECON.DM.setSegmentActs(ABER,a);
					%RECON.DM.setActs(ABER);
					
					F.planewave * RECON.A * RECON.DM;
                    if RECON.WFS.usePyr == 1
                        RECON.WFS.sensePyramid(F);
                    else
                        RECON.WFS.sense(F);
                    end
                    % RECON.WFS.quiver();
                    
                    % RECON.DM.show();
                    % F.show;
                    % hold on;
                    % RECON.WFS.quiver;
                    % hold off;
                    % bigtitle(sprintf('DH mode(%d,%d)',n,m));
                    % drawnow;
                    [x,y]=coords(F);
                    clf;
                    hold on;
                    subplot(1,2,1); 
                    quiver(RECON.WFS,1);
                    %pause;
                    subplot(1,2,2); 
                    imagesc(x,y,RECON.A.grid .* RECON.DM.grid); daspect([1 1 1]);
                    hold off;
                    drawnow;
					
					RECON.ACTS(:,nmode) = RECON.DM.actuators(:,3);
					RECON.SLOPES(:,nmode) = RECON.WFS.slopes;
					nmode = nmode + 1;
				end
				%subplot(1,2,1);imagesc(RECON.ACTS);bigtitle('Actuators');
				%subplot(1,2,2);imagesc(RECON.SLOPES);bigtitle('Slopes');
				%drawnow;
                % fprintf('\n');
			end
            % fprintf('\n');
			
            
            
            end
            
			
			fprintf('AOReconstructor: Using SVD to examine the slopes...\n');
			
			[UU,SS,VV] = svd(RECON.SLOPES','econ');
			RECON.U = UU;
			RECON.s = diag(SS);
			RECON.V = VV;
			
			clear UU SS VV
			
			RECON.rebuild;
		end
        
        % This tries to improve a reconstructor.  It doesn't matter how you
        % got the first one.
        function RECON = postProcess(RECON,D,Nmax,lambda)
			% RECON = zprogram(RECON,D,Nmax,lambda)
			
			if(nargin<5)
				if(isempty(RECON.lambda))
					RECON.lambda = AOField.RBAND;
				end
			else
				RECON.lambda = lambda;
            end
            
            % Save the original one.
            RECON0 = RECON.RECONSTRUCTOR;
            
            RECON.TrainingMethod = 'zernike postprocessed';
			RECON.OWD = Nmax;
			RECON.D = D;
			
			RECON.A.trueUp; 
			RECON.WFS.initBias(RECON.A);
			
			NZmodes = (Nmax+1)*(Nmax+2)/2;
			
			% We must do this for sin and cos.
			RECON.ACTS = zeros(RECON.DM.nActs,NZmodes);
			RECON.SLOPES = zeros(2*RECON.WFS.nSubAps,NZmodes);
			% 			disp(RECON);
			
			F = AOField(RECON.A);
			F.lambda = RECON.lambda;
			
			ABER = AOScreen(RECON.A);
			
			amp = RECON.amplitude;
			
			fprintf('AOReconstructor: Using ZERNIKES to feel out phase space to Order %d...\n',Nmax);
            fprintf('Number of modes being explored:%d\n',NZmodes);
			
			nmode = 1;
			
			for n=1:Nmax
                fprintf('\nOrder %d:',n);
                amp = 1*RECON.lambda/4/n;
				for m=-n:2:n
					ABER.zero.addZernike(n,m,amp,RECON.D);
					RECON.DM.setActs(ABER);
					F.planewave * RECON.A * RECON.DM;
                    %Added RECON.A 20091107. Was line below previously
					%F.planewave * RECON.DM;
                    if RECON.WFS.usePyr == 1
                        RECON.WFS.sensePyramid(F);
                    else
                        RECON.WFS.sense(F);
                    end
					RECON.ACTS(:,nmode) = RECON.DM.actuators(:,3);
					RECON.SLOPES(:,nmode) = RECON.WFS.slopes;
					nmode = nmode + 1;
                    fprintf('%d ',m);
                    [x,y]=coords(F);
%                     clf;
%                     hold on;
%                     subplot(1,2,1); 
%                     quiver(RECON.WFS,1);
                    %pause;
%                     subplot(1,2,2); 
%                     imagesc(x,y,RECON.A.grid .* RECON.DM.grid); daspect([1 1 1]);
%                     hold off;
				drawnow;
                    %pause;
				end
				%subplot(1,2,1);imagesc(RECON.ACTS);title('Actuators');
				%subplot(1,2,2);imagesc(RECON.SLOPES);title('Slopes');
                %drawnow;
                %fprintf('\n');
			end
			fprintf('\n');

            % patch up the experimental data.  The slopes need to be
            % prepended with the identity and the actuators with the
            % original reconstructor.
            
            RECON.ACTS = [RECON0 RECON.ACTS];
            RECON.SLOPES = [eye(RECON.WFS.nSubAps*2) RECON.SLOPES];
            
			fprintf('AOReconstructor: Using SVD to examine the slopes...\n');
            
            RECON.SLOPES(isnan(RECON.SLOPES)) = 0;
            RECON.ACTS(isnan(RECON.SLOPES)) = 0;
            
			[UU,SS,VV] = svd(RECON.SLOPES','econ');
			RECON.U = UU;
			RECON.s = diag(SS);
			RECON.V = VV;
			
			clear UU SS VV
			
			RECON.rebuild;
            
            clf;
            subplot(1,2,1);
            imagesc(RECON0);
            title('Original Reconstructor');
            colorbar;
            
            subplot(1,2,2);
            imagesc(RECON.RECONSTRUCTOR);
            title('New Reconstructor');
            colorbar;
            drawnow;
		end

        function WFS_SLOPES = processTestVectors(RECON,TESTVECTORS)
            % TESTVECTORS is a list of ACTUATOR positions.
            
            NVECTORS = size(TESTVECTORS,2);
            nSubAps = RECON.WFS.nSubAps;
            %nActs = RECON.DM.nActs;
            
            F = AOField(RECON.A);
			F.lambda = RECON.lambda;
            
            WFS_SLOPES = nan(2*nSubAps,NVECTORS);
            [Xwfs,Ywfs] = RECON.WFS.COORDS;
            BBox = RECON.A.BBox;
            BBox(1,:) = BBox(1,:) - 0.5;
            BBox(2,:) = BBox(2,:) + 0.5;
            
            for n=1:NVECTORS
                fprintf('%d ',n);
                RECON.DM.setActs(TESTVECTORS(:,n));
                F.planewave * RECON.A * RECON.DM;
                if RECON.WFS.usePyr == 1
                    RECON.WFS.sensePyramid(F);
                else
                    RECON.WFS.sense(F);
                end
                WFS_SLOPES(:,n) = RECON.WFS.slopes;
                
                if(RECON.verbose)
                    Xslopes = reshape(WFS_SLOPES(1:nSubAps,n),size(Xwfs));
                    Yslopes = reshape(WFS_SLOPES((1:nSubAps)+nSubAps,n),size(Xwfs));
                    RECON.A.show;
                    hold on;
                    quiver(Xwfs,Ywfs,Xslopes,Yslopes);axis xy;sqar;
                    xlim(BBox(:,1));% x and y may be reversed. :-0
                    ylim(BBox(:,2));
                    drawnow;
                    hold off;
                end
            end
            fprintf('\n');
        end
        
        function RECON = adhocProgram(RECON,width)
            
            RECON.ACTS = [];
            RECON.SLOPES = [];
            RECON.U = [];
            RECON.s = [];
            RECON.V = [];
            RECON.Nmodes = nan;
            RECON.TrainingMethod = 'jlc1';
            
            if(isempty(RECON.lambda))
                RECON.lambda = AOField.RBAND;
                fprintf('NOTE!!! I am setting the RECONSTRUCTOR lambda to %g.\n',...
                    RECON.lambda);
            end
            
            Nwfs = RECON.WFS.nSubAps;
            XY = RECON.WFS.subApCoords;
            xact = RECON.DM.actuators(:,1);
            yact = RECON.DM.actuators(:,2);
            
            amp = 1e-6; % TODO: Compute this!!!
            
            RECON.RECONSTRUCTOR = zeros(size(RECON.DM.actuators,1),2*Nwfs);
            ABER = AOScreen(RECON.DM);
            
            for direction=1:2
                for n=1:Nwfs
                    flag = 3 - direction; % sorry.
                    ABER.zero.addGDelta(XY(n,:),amp,width,flag);
                    RECON.RECONSTRUCTOR(:,n+(direction-1)*Nwfs) = ABER.interpGrid(xact,yact);
                    %RECON.show;drawnow;
                end
            end
                RECON.show;
        end
        
        function RECON = adhocProgram2(RECON,width)
            
            RECON.ACTS = [];
            RECON.SLOPES = [];
            RECON.U = [];
            RECON.s = [];
            RECON.V = [];
            RECON.Nmodes = nan;
            RECON.TrainingMethod = 'jlc1';
            
            if(isempty(RECON.lambda))
                RECON.lambda = AOField.RBAND;
                fprintf('NOTE!!! I am setting the RECONSTRUCTOR lambda to %g.\n',...
                    RECON.lambda);
            end
            
            Nwfs = RECON.WFS.nSubAps;
            XY = RECON.WFS.subApCoords;
            xact = RECON.DM.actuators(:,1);
            yact = RECON.DM.actuators(:,2);
            
            amp = 1e-6; % TODO: Compute this!!!
            
            RECON.RECONSTRUCTOR = zeros(size(RECON.DM.actuators,1),2*Nwfs);
            ABER = AOScreen(RECON.DM);
            
            for direction=1:2
                for n=1:Nwfs
                    flag = 3 - direction; % sorry.
                    ABER.zero.addGDelta(XY(n,:),amp,width,flag);
                    RECON.RECONSTRUCTOR(:,n+(direction-1)*Nwfs) = ABER.interpGrid(xact,yact);
                    %RECON.show;drawnow;
                end
            end
                RECON.show;
		end
		
		function RECON = rebuild(RECON,cutoff)
			if(nargin>1)
				RECON.Scutoff = cutoff;
            end
			
            if(isempty(RECON.s))
                fprintf('WARNING: RECONSTRUCTOR not trained using an SVD method.\n');
                %return;
                [UU,SS,VV] = svd(RECON.SLOPES','econ');
                RECON.U = UU;
                RECON.s = diag(SS);
                RECON.V = VV;
                
                clear UU SS VV
            end

			fprintf('The normalized singular value current cutoff is %g.\n', RECON.Scutoff);
			
			if(RECON.Scutoff<1)
				MODES = (RECON.s/RECON.s(1) >= RECON.Scutoff);
				RECON.Nmodes = sum(MODES);
			else
				if(RECON.Scutoff>length(RECON.s))
					RECON.Scutoff = length(RECON.s);
				end
				
				MODES = 1:round(RECON.Scutoff);
				RECON.Nmodes = MODES(end);
			end
			
			fprintf('This means that you will be correcting %d modes.\n',RECON.Nmodes);
			
			MODES(RECON.s(MODES)<1e-16) = [];
			
			RECON.RECONSTRUCTOR = ((RECON.V(:,MODES) * ...
				(diag(1./RECON.s(MODES)) * RECON.U(:,MODES)')) ...
				* RECON.ACTS')' ;

		end

		function FULL_RECON_MATRIX = fullRecon(RECON)
			NWFS = prod(size(RECON.WFS));
			NACT = RECON.DM.nActs;
			LIST = RECON.WFS.validSubapMap;
			LIST = [LIST;LIST+NWFS];
			
			FULL_RECON_MATRIX = zeros(NACT,2*NWFS);
			for n=1:NACT
				FULL_RECON_MATRIX(n,LIST) = RECON.RECONSTRUCTOR(n,:);
			end
		end
				
        function this = show(this)
            %figure;
            colormap(gray);
            imagesc(this.RECONSTRUCTOR);
            ylabel('actuators');
            xlabel('slopes');
            drawnow;
        end
        
        function wavenumber = k(this)
            wavenumber = 2*pi/this.lambda;
        end
    end
end

