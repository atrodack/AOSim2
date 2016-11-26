% An example of how to make a continuous face-sheet Deformable Mirror
%
% 20150412 ATRodack

clear all;
clc;
close all;

%% BMC DM
%Parameters
lambda = AOField.RBAND; % Red light.
SPACING = 1e-5; % fine spacing

%DM Specs
nActs = 1020; %32x32 minus 4 in the corners
Max_Stroke = 1.5e-6;
Pitch = 225e-6;
sidelength = (3.5)*10^-3;

%Coordinate System
xmin = -sidelength;
xmax = sidelength;
BMC_x = (xmin:SPACING:xmax);
ymin = -sidelength;
ymax = sidelength;
BMC_y = (ymin:SPACING:ymax);
[BMC_X,BMC_Y] = meshgrid(BMC_x,BMC_y);

BMC_pupil = ones(size(BMC_X));

%Construct the Pupil
Seg = AOSegment(length(BMC_pupil));
Seg.spacing(SPACING);
Seg.name = 'BMC Pupil';
Seg.grid(BMC_pupil);

%Make it an AOAperture Class
A = AOAperture;
A.spacing(SPACING);
A.name = 'BMC Aperture';
A.addSegment(Seg);
A.show;
drawnow;

%Make it an AODM
BMC_DM = AODM(A);
[X,Y] = BMC_DM.COORDS;
BMC_DM.show;
colormap(gray);
title('BMC Pupil');
drawnow;
pause(2);


%% Create Actuator Locations
%Make Actuator Coordiante Space
actuator_x = xmin:Pitch:xmax;
actuator_y = ymin:Pitch:ymax;
[ACTUATOR_X,ACTUATOR_Y] = meshgrid(actuator_x,actuator_y);
%Turn off Actuators at corners
ACTUATOR_X(1,1) = 0; ACTUATOR_X(32,32) = 0; ACTUATOR_X(32,1) = 0; ACTUATOR_X(1,32) = 0;
ACTUATOR_Y(1,1) = 0; ACTUATOR_Y(32,32) = 0; ACTUATOR_Y(32,1) = 0; ACTUATOR_Y(1,32) = 0;
ACTUATOR_X(ACTUATOR_X==0) = [];
ACTUATOR_Y(ACTUATOR_Y==0) = [];
%Write the locations so AOSim2 will understand them
BMC_ACTS = zeros(1020,2);
BMC_ACTS(:,1) = ACTUATOR_X(:);
BMC_ACTS(:,2) = ACTUATOR_Y(:);

% Add Actuators
BMC_DM.addActs(BMC_ACTS);
% Define Boundary Conditions
BMC_DM = BMC_DM.defineBC(sidelength+Pitch,5,'square');
% Plot Actuator Locations
BMC_DM.plotActuators;
title('BMC Actuator Locations');
pause(2);
% Plot Actuator Influence Regions
BMC_DM.plotRegions;
title('BMC Actuator Influence Regions');
pause(2);

%% Do something with the DM
% Flatten the Mirror
BMC_DM.flatten; %sets all actuator heights to zero


%% Make a WFS
d = sidelength;
WFS = AOWFS(A,d/9); %give it the aperture object copy properties, and the subaperture size
WFS.name = 'BMC Shack-Hartmann WFS';

%% Build the Reconstructor
RECON = AOReconstructor(A,BMC_DM,WFS); %inputs are the aperture, the DM, and the WFS
RECON.verbose = true; %turns on plotting

%% The Zernike Training Method (see zprogram method in AOReconstructor)
%This is exactly the zprogram method, just done out here to get a feel for
%what is happening to make WFS measurements and set DM positions. No AO is
%done. The Reconstructor Matrix is plotted at end at prompt from user

%Inputs to zprogram
Nmax = 9;%max Zernike Order
D = 1.05*d; %diameter of the Zernike pupil (make a little larger to avoid edge issues of Zernikes)

RECON.TrainingMethod = 'zernike';
RECON.OWD = Nmax; 
RECON.D = D; 

RECON.A.trueUp;

%Compute the initial Bias on the WFS (must be done before the WFS can be
%used, think of this as the planewave calibration to see where all your
%centroids are before adding crinkled light to the system)
RECON.WFS.initBias(RECON.A);

%Convert Zernike Order into # of Modes
NZmodes = (Nmax+1)*(Nmax+2)/2;

% We must do this for sin and cos. Sets up the matrices to store data in
RECON.ACTS = zeros(RECON.DM.nActs,NZmodes);
RECON.SLOPES = zeros(2*RECON.WFS.nSubAps,NZmodes);

% Make a WFS Field
F = AOField(RECON.A);
F.lambda = RECON.lambda;

% Make an aberration screen to set DM to
ABER = AOScreen(RECON.A);

amp = RECON.amplitude;

fprintf('AOReconstructor: Using ZERNIKES to feel out phase space to Order %d...\n',Nmax);
fprintf('Number of modes being explored:%d\n',NZmodes);

% Do the Programming over all modes
nmode = 1;
for n=1:Nmax
    fprintf('\nOrder %d:',n);
    amp = 1*RECON.lambda/4/n;
    for m=-n:2:n
        %Create the Zernike and set it to the AOScreen Object (the DM needs
        %it to be an AOScreen, or AOAtmo object to use setActs or setDM)
        weight = (n^2 + m^2).^(-5/6);
        ABER.zero.addZernike(n,m,weight*amp,RECON.D);
        % Set the DM face-sheet to the weighted Zernike
        RECON.DM.setActs(ABER);
        % Send a planewave through to get crinkled
        F.planewave * RECON.A * RECON.DM;
        
        % Make a WFS measurement (usePyr==1 uses a Pyramid Sensor Scheme,
        % else is a Shack-Hartmann. To set this, set the usePyr flag in the
        % WFS objected before putting into RECON, or after I suppose:
        % WFS.usePyr = 1; or WFS.usePyr = 0;
        if RECON.WFS.usePyr == 1
            RECON.WFS.sensePyramid(F);
        else
            RECON.WFS.sense(F); %Look in AOWFS.m to see what is happening and bask in its beauty
        end
        
        % Put the actuator positions and slopes into the initialized matrices
        RECON.ACTS(:,nmode) = RECON.DM.actuators(:,3);
        RECON.SLOPES(:,nmode) = RECON.WFS.slopes;
        nmode = nmode + 1;
        fprintf('%d ',m);
        [x,y]=coords(F);
        
        %Plot the Measured Slopes and Corresponding DM shape
        if(RECON.verbose)
            [x,y]=coords(F);
            clf;
            hold on;
            subplot(1,2,1);
            quiver(RECON.WFS,1);
            title(sprintf('Wavefront Sensor Slopes\n'));
            %pause;
            subplot(1,2,2);
            imagesc(x,y,RECON.A.grid .* RECON.DM.grid); daspect([1 1 1]);
%             RECON.DM.plotActuators;
            title(sprintf('DM Shape at mode %d\n',nmode));
            hold off;
            drawnow;
        end
    end
end
fprintf('\n');

fprintf('AOReconstructor: Using SVD to examine the slopes...\n');
% Clean up any NaNs
RECON.SLOPES(isnan(RECON.SLOPES)) = 0;
RECON.ACTS(isnan(RECON.ACTS)) = 0;
% Do SVD
[UU,SS,VV] = svd(RECON.SLOPES','econ');
RECON.U = UU;
RECON.s = diag(SS);
RECON.V = VV;

clear UU SS VV

% See the rebuild method in AOReconstructor
RECON.rebuild;

input('Press Enter to Continue');
%% Rebuild it with modes up to the max cutoff mode from SVD
RECON.rebuild(500);
% Show the Reconstructor Matrix and the log-scale plot of singular values
% used
subplot(1,2,1)
RECON.show;
title('Reconstructor Matrix');
subplot(1,2,2)
semilogy(RECON.s/RECON.s(1));
title(sprintf('Plot of the Singular Values\n up to the cutoff'))
display('Go look at figure 1 again!');

