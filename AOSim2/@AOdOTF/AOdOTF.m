classdef AOdOTF < AOField
    %AOdOTF A simple addition to AOSim2 that builds a system to use dOTF to
    %compute the wavefront phase.
    
    properties
        finger;
        PSF0;
        PSF1;
        OTF0;
        OTF1;
        dOTF;
        Phase;
        A;
        FOV;
        Plate_Scale;
        Mask;
        plotMask;
        Mask_interped;
        pupil_center;
        pupil_radius;
        calibration = true;
    end
    
    
    
    methods
        
        % Constructor
        function obj = AOdOTF(A,FOV,Plate_Scale)
            obj = obj@AOField(A);
            if nargin == 1
                obj.A = A;
                thld = (obj.lambda / (A.segList{1}.Segment.pupils(1,3)) * 206265);
                obj.FOV =  thld * 25; %make this kinda large to improve dOTF resolution
                obj.Plate_Scale = thld / 3;
            elseif nargin == 3
                obj.A = A;
                obj.FOV = FOV;
                obj.Plate_Scale = Plate_Scale;
            else
                error('Incorrect Inputs');
            end
                
        end
        
        function AOdOTF = calibrateWFS(AOdOTF,y_pos,x_pos,width,Field,ps)
            % This method calbirates the dOTF calculations.  It will create
            % a finger to place in the pupil, and then use a planewave on
            % the pupil to generate the dOTF.  The magnitude of this dOTF
            % is plotted, and used to create the masks needed for phase
            % retrieval. 3 points are selected, the center of the upper
            % pupil, the edge of the upper pupil, and the center of the
            % lower pupil.  The rest is done automatically.  Finally, all
            % of the dOTF calculating properties are cleared.
            
            if nargin == 1
                fprintf('Using Default Values\n');
                [x,y] = AOdOTF.A.coords;
                x_pos = max(x)/2;
                y_pos = max(y)/2;
                width = 0.25;
                Field = AOField(AOdOTF.A);
                ps = 1;
            elseif nargin == 4
                fprintf('Using Default Field\n');
                Field =AOField(AOdOTF.A);
                ps = 1;
            elseif nargin == 5
                ps = 1;
            end
            AOdOTF.create_finger(y_pos,x_pos,width);
            AOdOTF.mkPSF(Field,ps);
            AOdOTF.mkOTF(Field,ps);
            AOdOTF.mkdOTF;
            AOdOTF.plotdOTFframe;
            AOdOTF.mkMask;
            AOdOTF.truephase;
            AOdOTF.calibration = false;
            AOdOTF.cleardOTF;
        end %calibrateWFS
        
        function AOdOTF = sense(AOdOTF,Field,ps)
            AOdOTF.mkPSF(Field,ps);
            AOdOTF.mkOTF(Field,ps);
            AOdOTF.mkdOTF;
            AOdOTF.truephase;
        end %sense
            
        
        function AOdOTF = create_finger(AOdOTF,y_pos,x_pos,width)
            %Creates a pupil finger. y_pos is the height of the finger,
            %x_pos is the x position of the pupil, and width is the width
            %of the finger.
            
            if nargin == 1
                [x,y] = AOdOTF.A.coords;
                x_pos = max(x)/2;
                y_pos = max(y)/2;
                width = 0.25;
            end
            Finger = AOSegment(AOdOTF.A);
            Finger.name = 'Finger';           
            [X,Y] = Finger.COORDS;
            Finger.grid(~(Y<-y_pos & abs(X-x_pos)<width));
            AOdOTF.finger = Finger;
        end %create_finger
        
        function AOdOTF = mkPSF(AOdOTF,Field,ps)
            %Makes PSFs and stores them in the properties.  The first uses
            %the unmodified pupil, the second uses the modified pupil.
            %This is only to be called when PSFs are not input by hand (as
            %in actual lab data images). ps must be AOATMO or AOScreen
            %class from AOSim2, and acts as the phase aberration over the
            %pupil. Field must be AOField class from AOSim2.
            
            if nargin == 2
                ps = 1;
            end
%                 if ~isa(ps,'AOAtmo')
%                     if ~isa(ps,'AOScreen')
%                         error('ps must be of type AOAtmo or AOScreen');
%                     end
%                 end
            if ~isa(Field,'AOField')
                error('Field must be of type AOField');
            end
            
            if isempty(AOdOTF.finger)
                AOdOTF.create_finger;
            end
            Field.FFTSize = 1024;
            Field * ps * AOdOTF.A;
            AOdOTF.PSF0 = Field.mkPSF(AOdOTF.FOV,AOdOTF.Plate_Scale);
            Field.touch;
            
            Field * AOdOTF.finger;
            AOdOTF.PSF1 = Field.mkPSF(AOdOTF.FOV,AOdOTF.Plate_Scale);
            Field.touch;
        end
                
        
        function AOdOTF = mkOTF(AOdOTF,Field,ps)
            %Makes OTFs and stores them in properties.  If called with no
            %arguments, it simply Fourier Transforms what is stored in the 
            %PSF properties.  If called with one or two arguments, and the
            %PSF properties are empty, it will use Field and/or ps to
            %create and store PSFs, and use those to compute the OTF
            
            if nargin == 1
                ps = 1;
                Field = 1;
            elseif nargin == 2
                if ~isa(Field,'AOField')
                    error('Field must be of type AOField');
                end
                ps = 1;
                if isempty(AOdOTF.finger)
                    AOdOTF.create_finger;
                end
                
                if isempty(AOdOTF.PSF0)
                    Field.FFTSize = 1024;
                    Field.planewave * ps * AOdOTF.A;
                    AOdOTF.PSF0 = Field.mkPSF(AOdOTF.FOV,AOdOTF.Plate_Scale);
                    Field.touch;
                end
                
                if isempty(AOdOTF.PSF1)
                    Field * AOdOTF.finger;
                    AOdOTF.PSF1 = Field.mkPSF(AOdOTF.FOV,AOdOTF.Plate_Scale);
                    Field.touch;
                end
            elseif nargin == 3
%                 if ~isa(ps,'AOAtmo')  
%                     if ~isa(ps,'AOScreen')
%                         error('ps must be of type AOAtmo or AOScreen');
%                     end
%                 end
                if ~isa(Field,'AOField')
                    error('Field must be of type AOField');
                end
                if isempty(AOdOTF.finger)
                    AOdOTF.create_finger;
                end
                
                if isempty(AOdOTF.PSF0)
                    Field.FFTSize = 1024;
                    Field.planewave * ps * AOdOTF.A;
                    AOdOTF.PSF0 = Field.mkPSF(AOdOTF.FOV,AOdOTF.Plate_Scale);
                    Field.touch;
                end
                
                if isempty(AOdOTF.PSF1)
                    Field * AOdOTF.finger;
                    AOdOTF.PSF1 = Field.mkPSF(AOdOTF.FOV,AOdOTF.Plate_Scale);
                    Field.touch;
                end
            end

            AOdOTF.OTF0 = ifftshift(fft2(fftshift(AOdOTF.PSF0)));
            AOdOTF.OTF1 = ifftshift(fft2(fftshift(AOdOTF.PSF1)));     
        end %mkOTF
        
        function AOdOTF = mkdOTF(AOdOTF)
            %Computes the dOTF and stores it as a complex numbered matrix.
            %It will also store the phase of the dOTF in the Phase
            %property.  If running on 64-bit Ubuntu, the unwrapphase method
            %can be uncommented, and the phase stored will be unwrapped
            %using the unwt method (See unwrapphase method)
            
            AOdOTF.dOTF = (AOdOTF.OTF0) - (AOdOTF.OTF1);
            AOdOTF.Phase = angle(AOdOTF.dOTF);
%             AOdOTF.unwrapphase;
        end %mkdOTF
        
        function plotdOTFframe(AOdOTF)
            %Plots the magnitude of the dOTF. Most useful in the
            %calibration of the masks, but can be called if desired.
            
            imagesc(abs(AOdOTF.dOTF));
            bigtitle('dOTF Magnitude',12);
            colormap(gray);
            daspect([1 1 1]);
            axis xy;
            axis off;
        end
        
        
        function AOdOTF = unwrapphase(AOdOTF)
            %Unwrapps Phase.  See function uwrap in AOSim2\utils for more
            %information.  Can only be run on Ubuntu 64-bit OS as compiled.
            AOdOTF.Phase = uwrap(AOdOTF.Phase,'unwt');
        end %unwrapphase
        
        function AOdOTF = truephase(AOdOTF)
            %Retrieves the part of the dOTF phase that is important. This
            %uses the mask to grab the upper pupil conjugate, shifts it to
            %the central pixel, clips the matrix to remove surrounding
            %zeros, and then interpolates it up to the size of the DM
            %pupil (assuming AOdOTF.A is the size of the DM).  This might
            %only work for fingers that cause the dOTF to conjugate so
            %there is an upper left pupil and a lower right pupil. Other
            %configurations are UNTESTED.
            
            phase = AOdOTF.Phase;
            mask = AOdOTF.Mask;
            shift_point = AOdOTF.pupil_center;
            radius = round(AOdOTF.pupil_radius);
            
            [sizex,sizey] = size(phase);
            center = [round(sizex/2),round(sizey/2)];
            shift = [shift_point(1) - center(2),shift_point(2) - center(1)];
            
            phase_ref = phase(center(1),center(2));
            AOdOTF.unwrapphase;
            phase = AOdOTF.Phase - phase_ref;
            
            
            AOdOTF.Phase = phase .* mask;
%             AOdOTF.unwrapphase;
%             phase = AOdOTF.Phase - phase_ref;
            
            AOdOTF.Phase = phase;
            
            phase = circshift(phase,-shift);
            
            phase = phase - phase_ref;
            
            
            % resize to edge of Pupil
            phase = phase(center(1)-radius-0:center(1)+radius+0,center(2)-radius-0:center(2)+radius+0);
            if AOdOTF.calibration == true
                mask = circshift(mask,-shift);
                mask = mask(center(1)-radius - 0:center(1)+radius+0,center(2)-radius-0:center(2)+radius+0);
                AOdOTF.Mask_interped = mask;

            end
            AOdOTF.Phase = phase;
            AOdOTF.resize_phase_to_Pupil;
            AOdOTF.Phase = AOdOTF.Phase .* AOdOTF.Mask_interped;
            


        end
        
        function AOdOTF = resize_phase_to_Pupil(AOdOTF)
            %Interpolates the stored phase to the size of the DM pupil
            mask = double(AOdOTF.Mask_interped);
            phase = AOdOTF.Phase;
            [sizex,sizey] = size(phase);
            [xx,yy] = AOdOTF.coords;
            x = linspace(min(xx),max(xx),sizex);
            y = x;
            xq = linspace(min(x),max(x),length(xx));
            yq = xq;
            [X,Y] = meshgrid(x,y);
            [Xq,Yq] = meshgrid(xq,yq);
            AOdOTF.Phase = interp2(X,Y,phase,Xq,Yq);
            if AOdOTF.calibration == true
                AOdOTF.Mask_interped = interp2(X,Y,mask,Xq,Yq);
            end
        end

        
        function points = calibrate_pupil_mask(AOdOTF)
            %A quick method used to grab and store points off of an image.  
            %It uses the points to calculate the radius of the dOTF pupil 
            %and stores this.
            
            fprintf('\n***Please select point at center of the Pupil***\n');
            central_point = pickPoint;
            fprintf('\n***Please select point at edge of the Pupil***\n');
            edge_point = pickPoint;
            fprintf('\n***Please select point at center of the other Pupil***\n\n');
            central_point2 = pickPoint;
            points = cell(1,3);
            points{1} = central_point;
            points{2} = edge_point;
            points{3} = central_point2;
            AOdOTF.pupil_center = central_point;
            AOdOTF.pupil_radius = sqrt((points{2}(1) - points{1}(1))^2 + (points{2}(2) - points{1}(2))^2);
        end %calibrate_pupil_mask
        
        function AOdOTF = mkMask(AOdOTF)
            %Creates Boolean logic masks for the dOTF. plotMask is the
            %union of masks for both pupil conjugates to make plotting the
            %unprocessed dOTF phase look nicer.  Mask is the important
            %calculation, and is the piece used to retrieve the wanted
            %phase information.
            
            points = AOdOTF.calibrate_pupil_mask;
            [sizex,sizey] = size(AOdOTF.OTF0);
            xx = linspace(1,sizex,sizex);
            yy = linspace(1,sizey,sizey);
            [X,Y] = meshgrid(xx,yy);
            
            RA = sqrt((X-points{1}(2)).^2 + (Y-points{1}(1)).^2);
            A = RA<=AOdOTF.pupil_radius;
            
            RB = sqrt((X-points{3}(2)).^2 + (Y-points{3}(1)).^2);
            B = RB<=AOdOTF.pupil_radius;
            
            mask = A+B;
            mask(mask>0)=1;
            AOdOTF.plotMask = mask;
            AOdOTF.Mask = A & ~B;  
        end %mkMask
        
        
        function AOdOTF = cleardOTF(AOdOTF)
            %clears some of the properties
            
            AOdOTF.PSF0 = [];
            AOdOTF.PSF1 = [];
            AOdOTF.OTF0 = [];
            AOdOTF.OTF1 = [];
            AOdOTF.dOTF = [];
        end %cleardOTF
        
        
        function AOdOTF = show(AOdOTF)
            % A simple plotting tool that plots the stored PSFs, OTFs, the
            % magnitude of the dOTF, and the phase of the dOTF
            
            subplot(2,2,1)
            psf0 = AOdOTF.PSF0;
            psf1 = AOdOTF.PSF1;
            imagesc([(psf0).^0.5,(psf1).^0.5]);
            axis xy;
            axis off;
            bigtitle('PSF',12);
            daspect([1,1,1]);
            
            subplot(2,2,2)
            otf0 = AOdOTF.OTF0;
            otf1 = AOdOTF.OTF1;
            imagesc([abs((otf0).^0.5),abs((otf1).^0.5)]);
            axis xy;
            axis off;
            bigtitle('OTF',12);
            daspect([1,1,1]);
            
            subplot(2,2,3)
            imagesc(abs(AOdOTF.dOTF) .* AOdOTF.plotMask);
            bigtitle('dOTF Magnitude',12);
            daspect([1 1 1]);
            axis xy;
            axis off;

            subplot(2,2,4)
            imagesc((AOdOTF.Phase .* AOdOTF.A.combined.grid_));
            bigtitle('dOTF Phase',12);
            daspect([1 1 1]);
            axis xy;
            axis off;
            colormap(gray)
        end %show
        
    end
    
end
