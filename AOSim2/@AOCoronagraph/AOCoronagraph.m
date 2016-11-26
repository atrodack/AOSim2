classdef AOCoronagraph < AOSegment
    % AOCoronagraph class.
    % Works basically like an AOAperture.
    % Multiplying by an AOCoronagraph carries you to the Lyot plane.
    % 20151224: JLCodona.  UA.SO.CAAO.AOSim2
    
    properties
        REMAPPED    = []; % For use with PIAA.  If null, no PIAA.
        PPMASK        = []; % Pupil Mask for reducing matrix rank.
        APODIZER    = [];
        FPM         = []; % Focal plane Phase mask
        FPTM        = []; % Focal plane Transmission mask
        APCMLC      = []; % Defines FP region that defines apodization.
        
        LYOT        = []; % If null, then use the AOCoronagraph grid.
        Fnumber     = 1;
        D           = [];
        FL          = [];  % Focal length
       
        CENTROID = [];
        
        verbose     = false; % print debugging info.

        % I may want to protect this later.
    
        % This propagates the pupil field from the masked pixels to the
        % FPM.  
        
        % The Fourier matrix that carries Fpp_ to Ffp.  
        % I plan on having it use the full extent of Fpp_ and Ffp.  
        % Use the masks to reduce the domain and range.  
        % The pupil mask is imprinted on Fpp_ so it doesn't have to be in
        % the Fmatrix.
        Fmatrix = []; 
        
        % I don't think I need this anymore.
        invFmatrix = [];
        s = [];
        antiLyot_U = []; 
        
        % switches
        DO_APODIZER = false;
        
        % AOField caches
        
        Fpp = [];  % reference to the input pupil plane field.
        Fpp_ = []; % This is the field after apodization. PIAA remaps the pupil so this may be morphed.
        
        Ffp = [];  % Reference to the field in the focal plane, before the FPM.
        Ffp_ = []; % Field after the FPM.  Now probably useless.

        Flyot = []; % The field at the Lyot plane, but before any Lyot Stop is applied.
        
        Fscience = []; % The field at the final focal plane.

        % SELECTION MASKS
        
        FPM_ASSIGNED = [];
                
    end
    
    methods
        
        % Constructor
        function CORO = AOCoronagraph(varargin)
            CORO = CORO@AOSegment(varargin);
            CORO.name = 'Coronagraph';
            
            CORO.initAPODIZER();
            CORO.setupPPMASK(1e-5);
        end
          
        %% 
        function CORO = setupPPMASK(CORO,thresh)
            % CORO = setupPPMASK(CORO,thresh)
            % Use grid>thresh defines PUPIL boolean mask.
            
            CORO.PPMASK = logical(CORO.grid >= thresh);
        end

        function CORO = initAPODIZER(CORO)
            % CORO.initAPODIZER();
            CORO.APODIZER = AOSegment(CORO);
            CORO.APODIZER.name = 'Pupil Apodizer';
            CORO.APODIZER.constant(1);
        end
        
        function CORO = setupFPM(CORO,FocalLength,lambdaRef,pixel_ld,FPM_MAX_ld)
            % CORO = setupFPM(CORO,FocalLength,lambda,pixel_ld,FPM_MAX_ld)
            % FocalLength: The coronagraph focal length to the FPM.
            % Fnumber: The Fnumber coming into the focal plane.
            % lambdaRef: The reference wavelength for computing the plate
            % scale.
            % pixel_ld: How big is the pixel in lambda/D
            % FPM_MAX_ld: How many lambda/D is the radius of the FPM region.
           
            % Derive and estimate params
            CORO.lambdaRef = lambdaRef;
            CORO.FL = FocalLength;
            CORO.D = CORO.estimateD;
            CORO.Fnumber = CORO.FL/CORO.D;
            
            LAMBDA_D = CORO.Fnumber * lambdaRef;
            PIXEL = LAMBDA_D * pixel_ld;  % Fp plate scale.
            
            %% Make Focal Plane Masks
            N = 2*ceil(FPM_MAX_ld / pixel_ld)+1;
            CORO.FPM = AOScreen(N);
            CORO.FPM.spacing(PIXEL);
            
            CORO.FPM.zero;
            CORO.FPM.name = 'FPM OPD';
            
            CORO.FPTM = AOSegment(CORO.FPM);
            
            CORO.FPTM.constant(1.0);  % All pixels transmit.
            CORO.FPTM.name = 'CFPM';
            
            CORO.FPM_ASSIGNED = false(CORO.FPM.size); % Nothing is assigned.
            
            CORO.APCMLC = AOSegment(CORO.FPTM);
            CORO.APCMLC.name = 'APCMLC FP Region';
            AXIS = CORO.APCMLC.AXIS_PIXEL;
            CORO.APCMLC.zero.setPixel(AXIS(1),AXIS(2),1); % No apodization.
        end
        
        %%
        function CORO = initF(CORO,lambda)
        % CORO.initF(lambda);
        % Compute the operator that takes Fp_ from PP to FP before FPM.
        % NOTE: I use the FULL extent of Fp_ and Ffp for Fmatrix.
        % You may want to downselect rows and columns later with boolean
        % filters.
        
            fprintf('Initializing the Fourier operator...\n');
            tic;
            if(nargin<2)
                lambda = CORO.lambdaRef;
            end
            
            [Xp,Yp] = CORO.COORDS;
            [Xf,Yf] = CORO.FPM.COORDS;
            
            % Note that this isn't the proper measure factor.
            %CORO.Fmatrix = (prod(CORO.spacing)/CORO.numel) * exp((1i*2*pi/lambda/CORO.FL)*...
            %    ( Xf(:)*Xp(:)' + Yf(:)*Yp(:)') );
            CORO.Fmatrix = exp((1i*2*pi/lambda/CORO.FL)*( Xf(:)*Xp(:)'+Yf(:)*Yp(:)'));
            %toc
            
            if(~isempty(CORO.Fpp))
                FIELD = CORO.Fpp.copy;
            else
                FIELD = AOField(CORO);
                FIELD.lambda = CORO.lambdaRef;
                CORO.Fpp = FIELD;
            end
            
            %fprintf('Normalizing...\n');
            
            FIELD.zero.setPixel(FIELD.AXIS_PIXEL(1),FIELD.AXIS_PIXEL(2),1);
            NORM2 = norm(CORO.Fmatrix' * (CORO.Fmatrix * CORO.Fpp.vec));
            CORO.Fmatrix = CORO.Fmatrix/sqrt(NORM2);
            toc
            
        end
        
        function CORO = mkInvF(CORO,thresh)
            % CORO = mkInvF(CORO,thresh)
            % I don't think I need this anymore.
            
            fprintf('Building the pseudo-inverse...\n');
            if(nargin<2)
                thresh = 1e-8;
            end
            
            tic;
            [U,S,V] = svd(CORO.Fmatrix,'econ');
            CORO.s = diag(S);
            CORO.invFmatrix = pseudoInv_rebuild2(U,CORO.s,V,CORO.s(1)*thresh);
            CORO.antiLyot_U = U;
            toc
        end        
        
        %% Operators and tasks
        
        function CORO = PPtoFP(CORO,F)
            % CORO = CORO.PPtoFP([FIELD0])
            % Use Fmatrix to find Ffp starting from Fpp_.
            % Make sure Ffp_ is defined and properly initialized before
            % calling this.
            %
            % Note that I use literal element-by-element products of
            % same-sized grids to avoid interpolate errors.
            
            if(nargin==2)
                CORO.Fpp = F.copy;
                CORO.Fpp.name = 'Pupil Plane Field';
                CORO.Fpp.grid(F.grid);
            end
            
            CORO.Fpp_ = CORO.Fpp.copy;
            CORO.Fpp_.name = 'Field post-pupil and apodizer';
            if(isempty(CORO.APODIZER))
                CORO.Fpp_.grid(CORO.Fpp.grid .* CORO.grid);
            else
                CORO.Fpp_.grid(CORO.Fpp.grid .* CORO.grid .* CORO.APODIZER.grid);
            end

            if(isempty(CORO.Ffp))
                CORO.Ffp = AOField(CORO.FPTM);
                CORO.Ffp.lambda =CORO.lambdaRef;
                CORO.Ffp.name = 'Pre-FPM focal plane field';
            end
            
            if(isempty(CORO.Fmatrix))
                CORO.initF;
            end
            CORO.Ffp.grid(CORO.Fmatrix * CORO.Fpp_.grid_(:));
        end
        
        %         function CORO = PPtoFPmask(CORO,usePPMask)
        %             % CORO = CORO.PPtoFP([usePPMask])
        %             % Use Fmatrix to find Ffp starting from Fpp_.
        %             % Make sure Ffp_ is defined and properly initialized before
        %             % calling this.
        %             % The optional flag usePPMask (defaults true) determines use of
        %             % PPMASK.
        %
        %             if(nargin<2)
        %                 usePPMask = true;
        %             end
        %
        %             if(usePPMask)
        %                 CORO.Ffp.grid(...
        %                     reshape(CORO.Fmatrix(:,CORO.PPMASK(:)) * CORO.Fpp_.grid_(CORO.PPMASK(:)),...
        %                     CORO.FPM.size));
        %             else
        %                 CORO.Ffp.grid(reshape(CORO.Fmatrix * CORO.Fpp_.grid_(:),CORO.FPM.size));
        %             end
        %         end
        
        function CORO = PPtoLP(CORO,Fpp)
            % CORO.PPtoLP(Fpp)
            % 
            % This sets up the internal fields for whatever use.  
            % 
            % WARNING: This breaks the usual design model for AOSim2.
            % All the resulting fields are kept in the CORO object.
            
            if(nargin>1) % Fpp is provided
                CORO.Fpp = Fpp.copy;
            else
                if(isempty(CORO.Fpp))
                    CORO.Fpp = AOField(CORO);
                end
            end
            CORO.Ffp.lambda = CORO.lambdaRef;

            CORO.Fpp_ = CORO.Fpp.copy; % For post=aperture and apodizer.
            CORO.Fpp_.name = 'post-APOD Field';
            
            CORO.Ffp = AOField(CORO.FPM); % Pre-FPM.
            CORO.Ffp.lambda = CORO.lambdaRef;
            CORO.Ffp.name = 'Focal Plane Field';

            %CORO.Ffp_ = CORO.Ffp.copy; % For post-FPM.
            %CORO.Ffp_.name = 'post-FPM Field';

            % Fields are set up and initialized.
            % Now do the coronagraph calculations...
            
            CORO.Fpp_ * CORO; % pass through pupil.

            if(~isempty(CORO.APODIZER))
                CORO.Fpp_ * CORO.APODIZER; % Apodize.
            end

            if(isempty(CORO.Fmatrix))
                CORO.initF(CORO.lambdaRef);
                %CORO.initF(CORO.Fpp_.lambda).mkInvF(1e-8);
            end
            
            %% Compute the focal plane fields
            
            CORO.Ffp.grid(CORO.Fmatrix * CORO.Fpp_.grid_(:));
            
            CORO.FPtoLP;
        end
        
        function CORO = FPtoLP(CORO)
            % CORO.FPtoLP();
            % Propagate the pre-FPM field through to the Lyot Plane (pre-Lyot stop).
            % NOTE: This breaks the old design experiments.  They have been
            % removed from AOCronagraph.
            % NOTE: This only uses FPTM (i.e. not FPM) ATM.

            if(isempty(CORO.Ffp_))
                CORO.Ffp_ = CORO.Ffp.grid(CORO.Ffp.grid).copy; % This keeps the Ffp from being modified.
                CORO.Ffp_.name = 'post FPTM Field';
            end
            CORO.Ffp_.grid(CORO.Ffp.grid); % This keeps the Ffp from being modified.
            %CORO.Ffp_ * CORO.FPM * CORO.FPTM; % Uses interp2???
            CORO.Ffp_.grid(CORO.Ffp_.grid .* CORO.FPTM.grid);  % WARNING: ignoring FPM 
            
            % Compute the field in the Lyot plane.
            
            CORO.Flyot = CORO.Fpp.copy;  % Start with the relayed field.
            CORO.Flyot.name = 'Lyot Plane Field';

            CORO.Flyot.grid(CORO.Fmatrix'*CORO.Ffp_.grid_(:));
        end

        
        function CORO = findAPCMLC(CORO,count)
            % CORO.findAPCMLC(count);
            % Find the CORO.APODIZER that goes with the CORO.APCMLC spot.
            % Initialize the CORO.APODIZER before starting or it will
            % continue iterating for count times.
            
            %CORO.APODIZER.constant(1);
            if(isempty(CORO.Fpp))
                CORO.Fpp = AOField(CORO);
                CORO.Fpp.lambda = CORO.lambdaRef;
            end

            for n=1:count
                CORO.PPtoFP(CORO.Fpp.planewave);
                CORO.Ffp_ = CORO.Ffp.copy;
                CORO.Ffp_ * CORO.APCMLC;
                
                if(isempty(CORO.Flyot))
                    CORO.Flyot = CORO.Fpp.copy;
                    CORO.Flyot.name = 'Lyot Field';
                end
                
                CORO.Flyot.grid(CORO.Fmatrix'*CORO.Ffp_.grid_(:)).plotC(2);
                
                CORO.APODIZER.grid(abs(CORO.Flyot.grid));
                CORO.APODIZER.normalize;
                
                if(CORO.APCMLC.verbosity)                    
                    CORO.showCoro(4);
                    drawnow;
                end
            end
        end
        
        function CORO = showCoro(CORO,gamma)
            % CORO = showCoro(CORO,gamma)
            
            if(nargin<2)
                gamma = 4;
            end
            
            N1 = 3; N2 = 2;
            frame = 0;
            
            clf
            frame=frame+1;
            subplot(N1,N2,frame);
            CORO.APODIZER.show;
            title('Apodizer');
            axis off;
            
            frame=frame+1;
            subplot(N1,N2,frame);
            CORO.Fpp_.plotC(1);
            title('F_{pp}');
            axis off;
            
            frame=frame+1;
            subplot(N1,N2,frame);
            CORO.Ffp.plotC(gamma);
            title('F_{fp} (before masks)');
            axis off;
            
            frame=frame+1;
            subplot(N1,N2,frame);
            CORO.FPTM.plotC(gamma);
            colorbar off;
            title('CFPM');
            axis off;
            
            frame=frame+1;
            subplot(N1,N2,frame);
            if(~isempty(CORO.Ffp_))
                CORO.Ffp_.plotC(gamma);
            end
            title('F_{fp} (after masks)');
            axis off;
            
            frame=frame+1;
            subplot(N1,N2,frame);
            if(~isempty(CORO.Flyot))
                CORO.Flyot.plotC(2);
                title('F_{Lyot}');
            else
                title('F_{Lyot} (undefined)');
            end
            axis off;
            drawnow;
            
        end
        
        
        
        
        %% Coronagraph operations
        
        function CORO = response(CORO,ANGLE0)
            % CORO.response(ANGLE0)
            % Amplitude 1 response through the corograph.
            
            CORO.Fpp.planewave(1,ANGLE0);
            CORO.PPtoFP.FPtoLP;

%             CORO.Ffp_ = CORO.Ffp.copy;
%             CORO.Ffp_ * CORO.FPTM;
%             
%             CORO.Flyot.grid(CORO.Fmatrix' * CORO.Ffp_.grid_(:));
            
            if(isempty(CORO.LYOT))
                CORO.Flyot * CORO;
            else
                CORO.Flyot * CORO.LYOT;
            end
            
            CORO.Fscience = CORO.Ffp.copy;
            CORO.Fscience.name = 'Science Cam';
            CORO.Fscience.zero;
            psi = CORO.Fscience.grid;
            %psi(CORO.PPMASK(:)) = CORO.Fmatrix * CORO.Flyot.grid_(:);
            psi(:) = CORO.Fmatrix * CORO.Flyot.grid_(:);
            CORO.Fscience.grid(psi);
        end
        
        
        %%
        function LPvectors = LyotContribs(CORO,FPSELECT,JustPP)
         % LPvectors = CORO.LyotContribs([FPSELECT],[JustPP])
         % FPSELECT defaults to ALL PIXELS in FPM
         % JustPP defaults to false.
         % The LP field vectors (no masks) for each pixel in FPSELECT.
         %
         % NOTE: This does not include the effects of the FPM.
         
            if(nargin<2)
                FPSELECT = true(CORO.FPM.size);
            end
         
            if(nargin<3)
                JustPP = false;
            end
            
            psiF = CORO.Ffp.grid_(FPSELECT(:));
            % psiF = psiF .* CORO.FPTM.grid_(FPSELECT(:));
            psiF = psiF .* exp(1i*CORO.Ffp.k * CORO.FPM.grid_(FPSELECT(:)));

            if(JustPP)
                LPvectors = CORO.Fmatrix(FPSELECT(:),CORO.PPMASK(:))';
            else
                LPvectors = CORO.Fmatrix(FPSELECT(:),:)';
            end
                
            for n=1:size(LPvectors,2)
                LPvectors(:,n) = LPvectors(:,n) * psiF(n);
            end
         
        end
        
        function CORO = reset(CORO,FPSELECT)
            % CORO.reset([FPSELECT]);
            % Reset the FPM to the starting value.
            
            if(nargin>1)
                if(size(FPSELECT) == CORO.FPM.size)
                    CORO.FPM_ASSIGNED = FPSELECT;
                else
                    fprintf('ERROR: FPSELECT must be the same size as the FPM');
                    return;
                end
            end
            
            CORO.FPM.zero;
            CORO.FPTM.constant(1);
        end
        
        function [VMODES,s] = LyotModes(CORO,FIELD0,ThreshNum)
            % [MODES,s] = CORO.LyotModes(FIELD0,[ThreshNum]);
            % Compute the FPM modes that maximize power in the Lyot stop
            % pixels selected by CORO.PPMASK.
            % MODES are the selected svd V-modes and s is the full list of singular values.
            % ThreshNum determines how many modes to return.
            % no entry returns all of the modes.
            % ThreshNum > 1 returns that many modes.
            % ThreshNum < 1 returns those modes where s/s(1) > ThreshNum.
            
            CORO.PPtoFP(FIELD0);
            
            fprintf('Folding in the focal plane field...\n');
            tic;
            MATRIX = CORO.RVmerge(CORO.Fmatrix',CORO.Ffp.grid_(:));
            %MATRIX = CORO.Fmatrix' * sdiag(CORO.Ffp.grid_(:));
            %spdiags(CORO.Ffp.grid_(:),0,n,n)
            toc

            fprintf('SVD of the FP2LP operator...\n'); tic;
            %[~,S,VMODES] = svd(MATRIX(CORO.PPMASK(:),:),'econ');
            [~,S,VMODES] = svd(MATRIX(CORO.PPMASK(:),:));
            toc
            s = diag(S);
            clear S
            
            if(nargin < 3)
                return;
            end
            
            if(ThreshNum > 1) 
                Nmax = min(round(ThreshNum),size(VMODES,2));
                VMODES = VMODES(:,1:Nmax);
            else
                VMODES = VMODES(:,s/s(1)>ThreshNum);
            end
        end
        
        function [VMODES,s,SELECT,MATRIX] = FFPModes(CORO,FIELD0,ThreshNum)
            % [MODES,s,FFPMASK,CFPM2FFP_MATRIX] = CORO.FFPModes(FIELD0,[ThreshNum]);
            % Compute the FPM modes that maximize deliver power to the
            % final focal plane (FFP).
            %
            % MODES are the selected svd V-modes and s is the full list of singular values.
            % ThreshNum determines how many modes to return.
            % no entry returns all of the modes.
            % ThreshNum > 1 returns that many modes.
            % ThreshNum < 1 returns those modes where s/s(1) > ThreshNum.
            %
            % A non-Guru shouldn't need the FFPMASK and CFPM2FFP_MATRIX.
            % Save time and memory by calling as 
            % [MODES,s] = CORO.FFPModes(FIELD0);
            
            fprintf('Setting up the coronagraph...\n');
            %CORO.Fpp.grid(FIELD0.grid);
            CORO.PPtoFP(FIELD0);
            
            SELECT = (normalize(CORO.Ffp.mag2)>0.02);
            SELECT = SELECT(:);
            
            fprintf('Building the FP1-to-LYOT operator...\n'); tic;
            MATRIX = AOGrid.RVmerge(CORO.Fmatrix',CORO.Ffp.grid_);
            toc
            
            fprintf('Folding in the Lyot stop.\n');
            if(isempty(CORO.LYOT))
                fprintf('...using the PUPIL mask as the Lyot stop.\n')
                tic;
                MATRIX = AOGrid.LVmerge(MATRIX,CORO.grid_(:));
                toc
            else
                fprintf('...using CORO.LYOT as the Lyot stop.\n')
                tic;
                MATRIX = AOGrid.LVmerge(MATRIX,CORO.LYOT.grid_(:));
                toc
            end
            
            fprintf('Including the propagation to the final focal plane...\n');
            tic;
            MATRIX = CORO.Fmatrix(SELECT,:) * MATRIX;
            toc

            fprintf('Done building the operator.\n')

            fprintf('SVD...\n');
            tic;
            %[~,S,VMODES] = svd(MATRIX,'econ');
            [~,S,VMODES] = svd(MATRIX,'econ');
            %[~,S,VMODES] = svd(MATRIX,'econ');
            toc
            s = diag(S);
            clear S
            
            if(nargin < 3)
                return;
            end
            
            if(ThreshNum > 1) 
                Nmax = min(round(ThreshNum),size(VMODES,2));
                VMODES = VMODES(:,1:Nmax);
            else
                VMODES = VMODES(:,s/s(1)>ThreshNum);
            end
        end
        
        
        %% Utilities
        function CORO = setCentroid(CORO)
            % CORO = CORO.setCentroid();
            % Set the CENTROID property for determining the optical axis.
            
            CORO.CENTROID = CORO.centroid;
        end
    end
    
    %% Static methods (utils that I want to keep in this context)
    methods(Static=true)
        
        function CTRANS = BBtrans(K,H)
            % CTRANS = AOCoronagraph.BBtrans(FREQLIST,HEIGHTS);
            % FREQLIST is the list of lambda0/lambda that you want.
            % Heights is a list of NxM heights in units of lambda0.
            % The calculation will include N pixels with M contributions each.
            
            H = sum(H,2);
            K = K(:);
            
            CTRANS = 0;
            for n=1:length(H)
                CTRANS = CTRANS + exp(1i*K*H(n));
            end
            
            CTRANS = CTRANS/length(H);
        end
        
    end    
    
end

