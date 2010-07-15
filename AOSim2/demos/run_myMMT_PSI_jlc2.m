%% MMTSim.  
% Johanan L. Codona, Steward Observatory, University of Arizona.
% 20090426 JLCodona: First-light version.
% 20090601: PSI Sim batch run.

% % Load in my pre-built MMT model.
% load data/MMTAO_Model_JLC_20090426
% load data/MMTAO_Model_withBadActs_JLC_20090427

% WFS_FPS = 550;
WFS_FPS = 527;
EXPOSURE = 1/15;

NLOOKS = 1;

% gain=1; % gain>2 is asking for trouble!
gain = 0.3; % gain>2 is asking for trouble!
GAMMA = 2;  % This is the gamma correction for the PSF image.
DEX = 4;

SCIENCE_WAVELENGTH = AOField.MBAND;

PLATE_SCALE = 0.048574; % Clio arcsecs
FOV = 32*PLATE_SCALE;

%% Set up aberration model...
ABER = AOScreen(A);
ABER.radius = D/2;
ABER.zero;

if(false)
    % defocus Z(2,0).
    ABER.addZernike(2,0,SCIENCE_WAVELENGTH/16);
    
    % astigmatism Z(2,2).
    ABER.addZernike(2,-2,SCIENCE_WAVELENGTH/16);
    % ABER.addZernike(2,2,SCIENCE_WAVELENGTH/16);
    
    % coma Z(3,1).
    ABER.addZernike(3,-1,SCIENCE_WAVELENGTH/32);
    ABER.addZernike(3,1,SCIENCE_WAVELENGTH/24);
    
    % trefoil Z(3,3).
    ABER.addZernike(3,3, - SCIENCE_WAVELENGTH/32);
    ABER.addZernike(3,-3,SCIENCE_WAVELENGTH/16);
    
    % spherical aberration Z(4,0)
    ABER.addZernike(4,0, - SCIENCE_WAVELENGTH/16);
    
    % higher-order spherical aberration Z(2*n,0)
    ABER.addZernike(6,0, - SCIENCE_WAVELENGTH/24);
    ABER.addZernike(8,0, - SCIENCE_WAVELENGTH/32);
    ABER.addZernike(10,0, - SCIENCE_WAVELENGTH/48);
end

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
    
for nLook=1:NLOOKS

    %% Define the Atmosphere model and winds aloft.
    ATMO = AOAtmo(A);
    
    WFlow = AOScreen(1024,0.17,500e-9);
    WFlow.name = 'Lower altitude turbulence';
    WFhigh = AOScreen(2048,0.20,500e-9);
    WFhigh.name = 'High altitude turbulence';
    
    ATMO.addLayer(WFlow,1000);
    ATMO.addLayer(WFhigh,8000);
    
    ATMO.layers{1}.Wind = [3 1];
    ATMO.layers{2}.Wind = [1 -1]*20;
    
    r0 = ATMO.totalFriedScale;
    th_scat = AOField.VBAND/r0*206265;
    
    fprintf('The total r0 is %f cm.\n',100*ATMO.totalFriedScale);
    fprintf('The seeing is %.2f arcsecs.\n',th_scat);
    
    % Turning this off is like using dynamic refocus.
    ATMO.GEOMETRY = false;
    
    ATMO.BEACON = GUIDE_STAR; % Set this so ATMO knows how to compute the wavefront.
    
    %% Create the WFS and Science AOField objects...
    Fwfs = AOField(A);
    Fwfs.lambda = RECON.lambda;  % The Reconstructor was calibrated at a certain wavelength.
    
    F = AOField(A);
    F.lambda = SCIENCE_WAVELENGTH;  % Pick your favorite Science Wavelength at the top.
    F.FFTSize = 512*[1 1]; % This needs to be HUGE for the GMT.
    PSF = F.mkPSF(.2,.2/100); %TODO: This is a bug workaround.  FIXME!
    thld = F.lambda/D*206265;
    % HALO = F.mkHALO(6*thld,thld/4);
    HALO = F.mkPSF(FOV,PLATE_SCALE);

    %% Set your Primary and adaptive secondary initial conditions...
    A.trueUp;
    DM.setActs(0);
    % DM.addRippleActs(.3*[1 1],500e-9,0);
    % DM.addRippleActs(.5*[-1.2 .4],300e-9,pi/2);
    % for n=1:5
    %     DM.addRippleActs(.5*randn(1,2),600e-9*rand,2*pi*rand);
    % end
    
    % These touches should no longer be needed, but they don't hurt.
    touch(DM);
    touch(A);
        
    NT = WFS_FPS*3; % 3 second runs.
    TIMES = demean((1:NT)/WFS_FPS);
    
    SLOPES = zeros(2*WFS.nSubAps,length(TIMES));    
    HCUBE = zeros([size(HALO) length(TIMES)]);
    
    %% "The clock has started..."
    for n=1:NT
        t = TIMES(n);
        ATMO.time = t;
        fprintf('Look %d: time = %.04g\n',nLook,t);
        
        %% This is the guts of the AO closed-loop integrating servo....
        ATMO.BEACON = GUIDE_STAR;
        WFS.sense(Fwfs.planewave*ATMO*A*DM);
        SLOPES(:,n) = WFS.slopes;
        DM.bumpActs(-gain*RECON.RECONSTRUCTOR * WFS.slopes);
        DM.removeMean;
        % That was it!
        
        %% Meanwhile, in the Science Band...
        ATMO.BEACON = SCIENCE_OBJECT;
        F.planewave*ATMO*A*DM;

        F*ABER;  % Non-common path error.
        touch(F); % bug?
        PSF = F.mkPSF(FOV,PLATE_SCALE);
        PSF(isnan(PSF)) = 0;
        HCUBE(:,:,n) = PSF;
        % HCUBE(:,:,n) = F.mkHALO(FOV,PLATE_SCALE);
        
    end
    GOP = GOPdate;
    save(sprintf('noNCP_%s_Look%03d.mat',GOP,nLook),...
        'GOP','HCUBE','FOV','PLATE_SCALE','SCIENCE_WAVELENGTH',...
        'SLOPES','WFS_FPS','TIMES');
end
