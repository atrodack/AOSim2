dTHETA = 0.1;
RING = 1.5;
AR = 0.25;
dAZ = pi/32;
FOV = 2;
% r0 = 0.2;

NXPIX = 3;
NYPIX = 2;



% make_the_LBTI_AO_jlc
% make_the_MMT_AO_jlc



[x,y] = A.coords;

A1Center = [0 -7.2085];
A2Center = [0  7.2085];




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

FOV = 1;
dFOV = .01;

ActsL = 0;
ActsR = 0;

for nt=1:length(TIMES)
    
    ATMO.time = TIMES(nt);
    
    
    
    
    
    
    
end

