gain = 1.0; % This is the AO gain.  Set to 0 for no AO.

dTHETA = 0.1;

FOV = 1.5;
dFOV = .01;

%  r0 = 0.2;

NXPIX = 2;
NYPIX = 1;

%% Define the extended source as a set of point sources.
RING = 0.3; % radius of the circumstellar ring
OBJ2STAR = 1.0;
NSOURCE = 0;
PA = 60/360*2*pi;

AR = 0.25;
dAZ = 2*pi/NSOURCE;
AZS=0:dAZ:(2*pi-1e-3);

THETAS = [0 0];
BRIGHTS = 1.0;

if(NSOURCE>0)
    BRIGHT_PARTS = OBJ2STAR * BRIGHTS(1)/NSOURCE;
    
    ROT = [cos(PA) sin(PA);-sin(PA) cos(PA)];
    for AZ=AZS
        THETA = RING*(ROT*[cos(AZ); AR*sin(AZ)]);
        THETAS(end+1,:) = THETA;
        BRIGHTS(end+1) = BRIGHT_PARTS;
    end
end
AMPS = sqrt(BRIGHTS);

%% Continue with run setup...
make_the_LBTI_AO_jlc
% make_the_MMT_AO_jlc

[x,y] = A.coords;

A1Center = [0 -7.2085];
A2Center = [0  7.2085];

%% We need 2 AODMs for efficiency
DML = AODM(A);
DMR = AODM(A);
DML.name = 'LEFT DM672';
DMR.name = 'RIGHT DM672';
DML.addActs(LBT_Actuators,1,A);
DMR.addActs(LBT_Actuators,1,A);
DML.defineBC(8,16); % A circle of 16 null points at 8m radius.
DMR.defineBC(8,16); % A circle of 16 null points at 8m radius.

%% Define the working and combining fields

Fwfs = AOField(A); % start centered on the origin.
Fwfs.lambda = AOField.RBAND;

Fscience = AOField(A); % start centered on the origin.
Fscience.lambda = AOField.MBAND;

Fcombined = AOField([10 23]./A.dx); % start centered on the origin.
Fcombined.lambda = Fscience.lambda;
Fcombined.FFTSize = 2048;


%% Define the Atmosphere model and winds aloft.
ATMO = AOAtmo(A);

WFlow = AOScreen(1024,0.17,500e-9);
WFlow.name = 'Lower altitude turbulence';
WFhigh = AOScreen(2048,0.20,500e-9);
WFhigh.name = 'High altitude turbulence';

ATMO.addLayer(WFlow,1000);
ATMO.addLayer(WFhigh,8000);

ATMO.layers{1}.Wind = [5 0];
ATMO.layers{2}.Wind = [1 -1]*30;

% Turning this off is like using dynamic refocus.
ATMO.GEOMETRY = false;


%% Guide star selection
SODIUM_LAYER = 90e3;

LGS_BEACON0 = [0 0 1] * SODIUM_LAYER;
LGS_BEACON = [0 1/206265 1] * SODIUM_LAYER;  % Offset by 1 arcsec.
STAR = [0 0 1e10];
LEO_TARGET = [0 0 400e3];
GEOSYNC_TARGET = [0 0 42e6];

% NGS CASE
GUIDE_STAR = STAR; % pick one.
SCIENCE_OBJECT = STAR; % pick one.

% Na LGS CASE
% GUIDE_STAR = LGS_BEACON;
% SCIENCE_OBJECT = STAR;

% Looking at GeoSynchronous satellites using LGS
% GUIDE_STAR = LGS_BEACON;
% SCIENCE_OBJECT = GEOSYNC_TARGET;

ATMO.BEACON = GUIDE_STAR; % Set this so ATMO knows how to compute the wavefront.

%% Define the AO operation

dt = 1e-3;
TIMES = -0.5:dt:0.5;

EXPOSURE = 25;

[PSF,thx,thy] = Fcombined.mkPSF(FOV,dFOV);
CUBE = zeros([size(PSF) floor(length(TIMES)/EXPOSURE)]);

CCD = 0;

ActsL = 0; % Initial values for the actuators.
ActsR = 0;

[x,y] = A.coords;

CCD = 0;
nExp = 1;
for nt=1:length(TIMES)
    % for nt=1:5
    ATMO.time = TIMES(nt);
    
    frame = 0;
    clf;
    
    %% AO Part (now separated from the science part so we can loop over science sub-sources.)
    
    %% Left eye AO and Science... (try not to compute DM twice)
    Fwfs.Offset = ALCenter; % Move to the proper location to pick up the phase
    Fwfs.planewave; % make sure to apply the planewave here to get the baseline phase.
    Fwfs*ATMO; % This now gets the atmo phase at the right place.
    Fwfs.Offset = Center; % move to the common location for analysys...
    DML.setActs(ActsL); % Set the DM to the last LEFT settings.
    Fwfs*DML; % This makes the DM OPD map too.
    Fwfs*A;
    slopesL = WFS.sense(Fwfs).slopes;
    dActsL = RECON.RECONSTRUCTOR * slopesL; % save results for a while
    ActsL = ActsL - gain*dActsL; % LEFT AO is ready for next time.
    
%     frame=frame+1;subplot(NYPIX,NXPIX,1);
%     imagesc(x,y,Fscience.interferometer(1));
%     daspect([1 1 1]);
%     axis xy;
%     %Fwfs.show;
%     WFS.quiver;
%     setFoV(5);
%     title('Left WFS');
    
    %% Right eye AO...
    Fwfs.Offset = ARCenter; % Move to the proper location to pick up the phase
    Fwfs.planewave; % make sure to apply the planewave here to get the baseline phase.
    Fwfs*ATMO; % This now gets the atmo phase at the right place.
    Fwfs.Offset = Center; % move to the common location for analysys...
    DMR.setActs(ActsR); % set up the DM for the RIGHT side.
    Fwfs*DMR; % This makes the DM OPD map too.
    Fwfs*A;
    slopesR = WFS.sense(Fwfs).slopes;
    dActsR = RECON.RECONSTRUCTOR * slopesR; % save results for a while
    ActsR = ActsR - gain*dActsR; % RIGHT AO is ready for next time.
    
%     frame=frame+1;subplot(NYPIX,NXPIX,2);
%     imagesc(x,y,Fscience.interferometer(1));
%     daspect([1 1 1]);
%     axis xy;
%     %Fwfs.show;
%     WFS.quiver;
%     setFoV(5);
%     title('Right WFS');
    
    
    
    IMAGE = 0;
    
    %% SCIENCE LIGHT
    
    for n=1:length(BRIGHTS)
        
        THETA = THETAS(n,:);
        AMP = AMPS(n);
        
        %% LEFT Science.
        
        Fscience.Offset = ALCenter;
        Fscience.planewave(AMP,THETA)*ATMO; % make sure to call planewave at the right place!
        Fscience.Offset = Center; % move to the common location for manipulation.
        Fscience*A*DML;
        Fscience.Offset = ALCenter; % move it back for adding to the combiner
        
        Fcombined.zero;
        Fcombined + Fscience;
        
        % pictures...
%         frame=frame+1;
%         subplot(NYPIX,NXPIX,3);
        %     imagesc(x,y,Fscience.interferometer(1));
        %     daspect([1 1 1]);
        %     axis xy;
%         Fscience.show;
%         WFS.quiver(1); % the 1 overplots on the previous image.
%         title('Left Science');
        
        %% Right Science...
        
        Fscience.Offset = ARCenter;
        Fscience.planewave(AMP,THETA)*ATMO; % make sure to call planewave at the right place!
        Fscience.Offset = Center; % move to the common location for manipulation.
        Fscience*A*DMR;
        Fscience.Offset = ARCenter; % move it back for adding to the combiner
        Fcombined + Fscience;
        
        % pictures...
%         frame=frame+1;
%         subplot(NYPIX,NXPIX,4);
        %     imagesc(x,y,Fscience.interferometer(1));
        %     daspect([1 1 1]);
        %     axis xy;
%         Fscience.show;
%         WFS.quiver(1); % the 1 overplots on the previous image.
%         title('Left Science');

        %% Display what's happening...
%         frame=frame+1;
%         subplot(NYPIX,NXPIX,4);
        %     imagesc(x,y,Fscience.interferometer(1));
        %     daspect([1 1 1]);
        %     axis xy;
%         Fscience.show;
%         WFS.quiver(1);
%         title('Right Science');
        
        %% Beam combiner
        
%         frame=NXPIX+1;
%         subplot(NYPIX,NXPIX,frame+[0 1]);
%         frame=frame+1; % I used 2 slots.
%         Fcombined.show;
        
        
        %% PSF
        [PSF,thx,thy] = Fcombined.mkPSF(FOV,dFOV);
        IMAGE = IMAGE + PSF;
        
%         drawnow;
    end % n star loop... produces IMAGE
    
    CCD = CCD + IMAGE;
    
    %imagesc(thx,thy,log10(normalize(PSF)),[-3 0]);
    frame=frame+1; subplot(NYPIX,NXPIX,frame);
    imagesc(thx,thy,log10(normalize(IMAGE)),[-3 0]);
    daspect([1 1 1]);
    axis xy;
    title('Fizeau Image');
    
    frame=frame+1; subplot(NYPIX,NXPIX,frame);
    imagesc(thx,thy,log10(normalize(CCD)),[-3 0]);
    daspect([1 1 1]);
    axis xy;
    title(sprintf('Fizeau Exposure: t=%.4g',ATMO.time));
    drawnow;
    
    if(mod(nt,EXPOSURE)==0)
        CUBE(:,:,nExp) = CCD;
        HEADER = struct();
        HEADER.GOPDATE = GOPdate;
        HEADER.AO_WFS = dt;
        HEADER.EXPOSURE = EXPOSURE;
        HEADER.RING = RING;
        HEADER.OBJ = OBJ2STAR;
        HEADER.NSOURCES = NSOURCE;
        HEADER.PA = PA;
        HEADER.nExp = nExp;
        nExp = nExp + 1;
        CCD = 0;
        fits_write_image('LBT_FIZEAU_extended_withAO.fits',CUBE,HEADER);
    end
    
end % loop over time

