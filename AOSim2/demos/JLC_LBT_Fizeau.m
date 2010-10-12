gain = -1;

dTHETA = 0.1;
RING = 1.5;
AR = 0.25;
dAZ = pi/32;
FOV = 2;
% r0 = 0.2;

NXPIX = 4;
NYPIX = 3;

% make_the_LBTI_AO_jlc
% make_the_MMT_AO_jlc

[x,y] = A.coords;

A1Center = [0 -7.2085];
A2Center = [0  7.2085];

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

CCD = 0;

FOV = 2;
dFOV = .01;

ActsL = 0;
ActsR = 0;

[x,y] = A.coords;

for nt=1:length(TIMES)
    frame = 0;
    
    ATMO.time = TIMES(nt);
    
    Fcombined.zero;

    %% Left eye AO and Science... (try not to compute DM twice) 
    DM.setActs(ActsL);
    Fwfs.Offset = ALCenter; % Move to the proper location to pick up the phase
    Fwfs.planewave; % make sure to apply the planewave here to get the baseline phase.
    Fwfs*ATMO; % This now gets the atmo phase at the right place.
    Fwfs.Offset = Center; % move to the common location for analysys...
    DM.setActs(ActsL); % This makes the DM OPD map too.  
    Fwfs*DM; % This makes the DM OPD map too.  
    Fwfs.planewave*ATMO; % This now gets the atmo phase at the right place.
    Fwfs*A;
    dActsL = RECON.RECONSTRUCTOR * WFS.sense(Fwfs).slopes; % save results for a while
    
    frame=frame+1;subplot(NYPIX,NXPIX,frame);
    %     imagesc(x,y,Fscience.interferometer(1));
    %     daspect([1 1 1]);
    %     axis xy;
    Fwfs.show;
    title('Left WFS');
    
    % Now the LEFT science while we have the DM set.
    
    Fscience.Offset = ALCenter;
    Fscience.planewave*ATMO;
    Fscience.Offset = Center; % move to the common location for manipulation.
    Fscience*A*DM;
    Fscience.Offset = ALCenter; % move it back for adding to the combiner
    Fcombined + Fscience;

    % pictures...
    frame=frame+1;subplot(NYPIX,NXPIX,frame);
    %     imagesc(x,y,Fscience.interferometer(1));
    %     daspect([1 1 1]);
    %     axis xy;
    Fscience.show;
    title('Left Science');
    
    ActsL = ActsL - gain*dActsL; % ready for next time.
    
    %% Right eye AO and Science... (try not to compute DM twice)
    DM.setActs(ActsR);
    Fwfs.Offset = ARCenter; % Move to the proper location to pick up the phase
    Fwfs.planewave; % make sure to apply the planewave here to get the baseline phase.
    Fwfs*ATMO; % This now gets the atmo phase at the right place.
    Fwfs*DM.setActs(ActsR); % This makes the DM OPD map too.  
    Fwfs.planewave*ATMO; % This now gets the atmo phase at the right place.
    Fwfs.Offset = Center; % move to the common location for pupil, DM, and analysys...
    Fwfs*A;
    dActsR = RECON.RECONSTRUCTOR * WFS.sense(Fwfs).slopes; % save results for a while
    
    frame=frame+1;subplot(NYPIX,NXPIX,frame);
    %     imagesc(x,y,Fscience.interferometer(1));
    %     daspect([1 1 1]);
    %     axis xy;
    Fwfs.show;
    title('Right WFS'); 
    
    
    % Now the RIGHT science field while we have the DM set.
    
    Fscience.Offset = ARCenter;
    Fscience.planewave*ATMO;
    Fscience.Offset = Center; % move to the common location for manipulation.
    Fscience*A*DM;
    Fscience.Offset = ARCenter; % move it back for adding to the combiner
    Fcombined + Fscience;

    ActsR = ActsR - gain*dActsR; % ready for next time.
    
    %% Display what's happening...
    frame=frame+1;subplot(NYPIX,NXPIX,frame);
    %     imagesc(x,y,Fscience.interferometer(1));
    %     daspect([1 1 1]);
    %     axis xy;
    Fscience.show;
    title('Right Science');
    
    frame=frame+1;subplot(NYPIX,NXPIX,frame);
%     figure(2);
    Fcombined.show;
%     figure(1);
    frame=frame+1;subplot(NYPIX,NXPIX,frame);
    PSF = Fscience.mkPSF(FOV,dFOV);
    
    imagesc(log10(normalize(PSF)),[-3 0]); 
    daspect([1 1 1]);
    axis xy;
    
    drawnow;
    
end

