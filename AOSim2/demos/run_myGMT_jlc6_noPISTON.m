%% GMTSim2.  
% Johanan L. Codona, Steward Observatory, University of Arizona.
% 20090426 JLCodona: First-light version.

% % Load in my pre-built MMT model.
% load data/MMTAO_Model_JLC_20090426
% load data/MMTAO_Model_withBadActs_JLC_20090427
% load data/MMTAO_Model_working

MAGIC_PISTONS = true;

%WFS = WFS17; % local alias.

WFS_FPS = 1000;

% gain=1; % gain>2 is asking for trouble!
% gain = 0.3; % gain>2 is asking for trouble!
gain = 0.5; % gain>2 is asking for trouble!
pgain = 0.1; % pgain is the piston gain.

GAMMA = 2;  % This is the gamma correction for the PSF image.
SCIENCE_WAVELENGTH = AOField.MBAND;

FOV_START = 0.3;
FOV_AO_ON = 0.3;

AO_STARTTIME = 0.0;
ZOOM_STARTTIME = 0.01;
ZOOM_ENDTIME = AO_STARTTIME;

%% Define the Atmosphere model and winds aloft.

% Height (m)  Cn2.dh (m^{1/3})  Wind speed (m/s)  Wind direction (deg.)
% 100        1.28e-13           10                0
% 500        9.88e-15           10                10
% 1000       2.11e-14           15                15
% 2000       5.54e-14           18                20
% 4000       3.77e-14           20                35
% 8000       3.39e-14           50                45
% 16000      3.94e-14           30                30

% Height (m)  Cn2.dh (m^{1/3})  Wind speed (m/s)  Wind direction (deg.)

ATMO_PROFILE = [
100        1.28e-13           10                0
% 500        9.88e-15           10                10
% 1000       2.11e-14           15                15
2000       5.54e-14           18                20
4000       3.77e-14           20                35
% 8000       3.39e-14           50                45
16000      3.94e-14           30                30
];

ATMO_PROFILE(:,5) = 2.^ceil(log2(round((D+ATMO_PROFILE(:,3)*2)/A.dx)));

ATMO = AOAtmo(A);

for n=1:size(ATMO_PROFILE,1)
	WEATHER = AOScreen(ATMO_PROFILE(n,5),0.2,500e-9);
	WEATHER.setCn2(ATMO_PROFILE(n,2));
	
	WEATHER.name = sprintf('Layer %d: alt %g',n,ATMO_PROFILE(n,1));;
	
	ATMO.addLayer(WEATHER,ATMO_PROFILE(n,1));
	
	az = pi*ATMO_PROFILE(n,4)/180;
	ATMO.layers{n}.Wind = [cos(az) sin(az)]*ATMO_PROFILE(n,3);
end

% WFlow = AOScreen(1024,0.17,500e-9);
% WFlow.name = 'Lower altitude turbulence';
% WFhigh = AOScreen(2048,0.20,500e-9);
% WFhigh.name = 'High altitude turbulence';
% 
% ATMO.addLayer(WFlow,1000);
% ATMO.addLayer(WFhigh,8000);
% 
% ATMO.layers{1}.Wind = [3 1];
% ATMO.layers{2}.Wind = [1 -1]*20;

r0 = ATMO.totalFriedScale;
th_scat = AOField.VBAND/r0*206265;

fprintf('The total r0 is %f cm.\n',100*ATMO.totalFriedScale);
fprintf('The seeing is %.2f arcsecs.\n',th_scat);

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

%% Create the WFS and Science AOField objects...
Fwfs = AOField(A);
Fwfs.lambda = RECON.lambda;  % The Reconstructor was calibrated at a certain wavelength.

F = AOField(A);
F.lambda = SCIENCE_WAVELENGTH;  % Pick your favorite Science Wavelength at the top.
F.FFTSize = 2048*[1 1]; % This needs to be HUGE for the GMT.
PSF = F.mkPSF(.2,.2/100); %TODO: This is a bug workaround.  FIXME!

%% Set your Primary and adaptive secondary initial conditions...

A.lambdaRef = F.lambda; % for plotting only.

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

clf;
colormap(gray);  % Looks more official.

% This is to select the points that will be included in the phase
% histogram.
mask = (A.grid>0.5);

% This is the brightest pixel seen to date.
Ipeak = 0;  
HISTMAX = 2000;
[x,y]=coords(F);

TIMES = (1:1500)/WFS_FPS;

ALL_PISTONS = zeros(7,length(TIMES));

% Strehl plot setup
mask = (A.grid>0.5);
STREHL = zeros(size(TIMES));
maxStrehl = 0.3;
N1=2;N2=3;  % Selects display geometry.

CCD = 0;

profile on

%% "The clock has started..."
for n=1:length(TIMES)
    t = TIMES(n);
    ATMO.time = t-TIMES(1000);
	
    %% This is the guts of the AO closed-loop integrating servo....
	
    ATMO.BEACON = GUIDE_STAR;
    WFS.sense(Fwfs.planewave*ATMO*A*DM);
    
    if(t>AO_STARTTIME)  % Suffer with seeing limit for 50ms.
        DM.bumpActs(-gain*RECON.RECONSTRUCTOR * WFS.slopes);
        DM.removeMean;
		
		% Don't do pistons at the WFS band, use longer wavelengths.
% 		if(MAGIC_PISTONS)
% 			PISTONS = WFS.magicPistonSensor(Fwfs,A);
% 			A.bumpPistons(-pgain*PISTONS);
% 		end
	end
    % That was it!
    	
    %% Meanwhile, in the Science Band...
    ATMO.BEACON = SCIENCE_OBJECT;
	F.planewave*ATMO*A*DM;
	
	if(t>AO_STARTTIME)  % Suffer with seeing limit for 50ms.
		% This is NEW!
		if(MAGIC_PISTONS)
			PISTONS = WFS.magicPistonSensor(F,A);
			A.bumpPistons(-pgain*PISTONS);
			ALL_PISTONS(:,n) = PISTONS;
		end
	end
	
    %% Plot some interesting pictures...
    clf; % Don't screw around.  Just clear it.
    
	subplot(N1,N2,1);
    
	if(t < ZOOM_STARTTIME)
		FOV = FOV_START;
	elseif(t < ZOOM_ENDTIME)
		FOV = (FOV_AO_ON-FOV_START)*...
			(t-ZOOM_STARTTIME)/(ZOOM_ENDTIME-ZOOM_STARTTIME)...
			+FOV_START;  % In arcsecs.
	else
		    FOV = FOV_AO_ON;  % In arcsecs.
	end

	RNG = FOV * [-1 1];
	PSF = F.mkPSF(FOV,FOV/100);
	
	if(t>ZOOM_ENDTIME) 
		CCD = CCD + PSF;
		if(mod(n,50)==0)
			fits_write_image('/tmp/GMTSim2_exposure.fits',CCD/max(CCD(:)));
		end
	end
	
	Ipeak = max(Ipeak,max(PSF(:)));
	imagesc(RNG,RNG,(PSF/Ipeak).^(1/GAMMA));
	axis xy;
    daspect([1 1 1]);
    title(sprintf('PSF (\\lambda=%.2g microns, \\gamma=%g) t=%.3f',...
        F.lambda*1e6,GAMMA,t));
    xlabel('arcsecs');
    ylabel('arcsecs');
    
	subplot(N1,N2,2);
	A.show;
    colorbar off;
	WFS.quiver(1);
    %title('WFS Slopes (autoscaled)');
    title('WFS Slopes');
	axis xy;
	
    subplot(N1,N2,3);
	
	[x,y]=F.coords;
    imagesc(x,y,(F.interferometer(.75)),[0 3]);
    title('Science Band Interferometer');
	daspect([1 1 1]);
	axis xy;
	% DM.plotActuators;
	
    %%
    subplot(N1,N2,4);
	g = F.grid_;
    STREHL(n) = abs(mean(g(mask)))^2;
	maxStrehl = min(1,max(maxStrehl,STREHL(n)*1.2));
    plot(TIMES(1:n),STREHL(1:n),'k-');
    ylim([0 maxStrehl]);
    xlim([0 TIMES(end)]);
    %xlim([-0.25 0]+t);
    ylabel('Strehl');
    xlabel('t (secs)');
	title(sprintf('Strehl History: Correcting %d modes, gain=%.2f.',...
		RECON.Nmodes,gain));
	clear g

	subplot(N1,N2,5:6);
	surf(x,y,DM.grid,abs(A.grid),'LineStyle','none');
% 	surf(x,y,DM.grid,'LineStyle','none');
	zlim([-1 1]*3e-6);
	daspect([1 1 5e-6]);
	lt=light();
	set(lt,'Position',[-1 0 1]);
    % You may need this if you aren't saving the frames.
    drawnow;

    %% This saves the current picture as a JPEG.
    filename = sprintf('/tmp/FRAME_%04d.jpg',n);
    rez = 160;

    resolution = sprintf('-r%d',rez);
    print(resolution,'-djpeg',filename);
	
	if(ATMO.time>2.0)
		break;
	end
end

