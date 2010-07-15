%% Johanan L. Codona's 10" LX200 Simulator. 
% (Note: This telescope does not have adaptive optics.)  
% Johanan L. Codona, Steward Observatory, University of Arizona.
% 20090827 JLCodona: First-light version.

WFS_FPS = 1000;
GAMMA = 2;  % This is the gamma correction for the PSF image.
SCIENCE_WAVELENGTH = AOField.VBAND;
% SCIENCE_WAVELENGTH = 400e-9;

D = 21*0.0254;
lambda0 = 500e-9;

FOV = 2.5;
dFOV = lambda0/D*206265/2;

FFT_SIZE = 256;

FPS = 30;
NEXP = round(WFS_FPS/FPS);

%% Define the Atmosphere model and winds aloft.
ATMO = AOAtmo(A);

WFlow = AOScreen(1024,0.10,500e-9);
WFlow.name = 'Lower altitude turbulence';
WFhigh = AOScreen(2048,0.15,500e-9);
WFhigh.name = 'High altitude turbulence';

ATMO.addLayer(WFlow,1000);
ATMO.addLayer(WFhigh,8000);

ATMO.layers{1}.Wind = [3 1];
ATMO.layers{2}.Wind = [1 -1]*20;

r0 = ATMO.totalFriedScale;
th_scat = AOField.VBAND/r0*206265;

fprintf('The total r0 is %.1f cm.\n',100*ATMO.totalFriedScale);
fprintf('The V-band seeing is %.2f arcsecs.\n',th_scat);

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

ATMO.BEACON = GUIDE_STAR; % Set this so ATMO knows how to compute the wavefront.

%% Create Science AOField objects...
F = AOField(A);
F.lambda = SCIENCE_WAVELENGTH;  % Pick your favorite Science Wavelength at the top.
F.FFTSize = FFT_SIZE*[1 1]; % This should be thought about for best performance.
PSF = F.mkPSF(FOV,dFOV); %TODO: This is a bug workaround.  FIXME!

%% Set your Primary and adaptive secondary initial conditions...
A.trueUp;
touch(A);

% This is to select the points that will be included in the phase
% histogram.
mask = (A.grid>0.5);

% This is the brightest pixel seen to date.
Ipeak = 0;  
HISTMAX = 2000;
[x,y]=coords(F);

N1=1;N2=2;  % Selects display geometry.

CCD = 0;
nCCD = 0;
nCCD_FRAME = 1;
CUBE = [];

ATMO.BEACON = SCIENCE_OBJECT;

%% "The clock has started..."
for n=1:2000
	% 1kHz Frame rate.  Divide by 2000 to match CoDR.
    t = n/WFS_FPS;
    ATMO.time = t-1.5;
	
    fprintf('%d:%d ',n,nCCD);
    
	F.planewave*ATMO*A;
	
	PSF = F.mkPSF(FOV,dFOV);
    
    CCD = CCD + PSF;
    nCCD = nCCD + 1;
    if(nCCD>=NEXP)
        nCCD = 0;
        CUBE(:,:,end+1) = CCD;
        CCD = 0;
        HEADER = struct;
        HEADER.sim = 'JLC 10in LX200 Sim';
        HEADER.lambda = SCIENCE_WAVELENGTH;
        HEADER.WFS_FPS = WFS_FPS;
        HEADER.FPS = FPS;
        HEADER.FOV = 2*FOV;
        HEADER.pixel = dFOV;
        r0 = ATMO.totalFriedScale;
        HEADER.r0 = r0;
        HEADER.seeing = AOField.VBAND/r0*206265;

        fits_write_image('/tmp/LX200_VideoSim.fits',CUBE,HEADER);
        fprintf('\n');

    end
end

